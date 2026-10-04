# -*- coding: utf-8 -*-
"""
thrombosis_data_fetcher.py
==========================
Сборщик данных по тромбозам из открытых источников + Физиологические профили.

Источники (в порядке приоритета):
1. WHO Mortality DB (CSV через прямую ссылку)
2. Kaggle datasets (через kaggle API, если установлен)
3. Встроенные базовые данные из клинических гайдлайнов
4. Гемодинамические профили для 0D CFD симуляции
"""

import os
import sys
from pathlib import Path

import pandas as pd
import requests

SOURCE = "all"
_in_jupyter = "ipykernel" in sys.modules
if not _in_jupyter and len(sys.argv) > 1:
    _valid = {"all", "who", "kaggle", "pubmed", "local"}
    SOURCE = sys.argv[1] if sys.argv[1] in _valid else SOURCE

OUT_DIR = Path("/mnt/c/Work/Cursor_CFD/data_science_aorta") if sys.platform != "win32" else Path("c:/Work/Cursor_CFD/data_science_aorta")
os.makedirs(OUT_DIR, exist_ok=True)

CSV_OUT      = OUT_DIR / "thrombosis_stats.csv"
ENRICHED_OUT = OUT_DIR / "thrombosis_stats_enriched.csv"
WHO_OUT      = OUT_DIR / "who_cvd_mortality.csv"
UCI_OUT      = OUT_DIR / "uci_heart_disease.csv"
BIOPHY_OUT   = OUT_DIR / "cfd_biophysical_triggers.csv"
RISK_OUT     = OUT_DIR / "who_cvd_risk_factors.csv"
HEMO_OUT     = OUT_DIR / "hemodynamic_profiles.csv"

# ==============================================================================
# 1. ЦЕЛЕВЫЕ КЛИНИЧЕСКИЕ ПОКАЗАТЕЛИ И МАППИНГ ПРИЧИН СМЕРТИ
# ==============================================================================

# Слой маппинга: Причина из CSV ВОЗ -> Клинический сценарий модели
CAUSE_TO_SCENARIO = {
    "Ischemic stroke (cardioembolic)": "stroke",
    "Carotid artery thrombosis":       "stroke",
    "Cerebral venous sinus thrombosis":"stroke",

    "Coronary thrombosis (STEMI)":     "heart_attack",
    "Atrial thrombus (LAA)":           "heart_attack",

    "Pulmonary embolism":              "cardiac_arrest",  # Массивная ТЭЛА ведет к остановке
    "Mesenteric artery thrombosis":    "cardiac_arrest",  # Шок и остановка
    "Aortic thrombus":                 "combined",        # Тяжелая патология аорты

    "Deep vein thrombosis":            "baseline",        # Само по себе сердце не убивает
    "Renal artery thrombosis":         "hypertension",    # Вызывает почечную гипертензию
}

CLINICAL_TARGETS = [
    # ХРОНИЧЕСКИЕ ПАТОЛОГИИ
    {
        "scenario": "baseline", 
        "HR": 72, "SBP": 120, "DBP": 80, "PWV": 7.0, 
        "compliance_change": 0.0, "pulse_pressure": 40, "CO": 5.0, "EF": 60,
        "IMT": 0.65, "AVA": 3.5, "peak_velocity": 1.2, "mean_gradient": 4.0,
        "clinical_notes": "Healthy normal"
    },
    {
        "scenario": "hypertension", 
        "HR": 76, "SBP": 140, "DBP": 90, "PWV": 9.5, 
        "compliance_change": -0.60, "pulse_pressure": 50, "CO": 5.0, "EF": 55,
        "IMT": 0.75, "AVA": 3.5, "peak_velocity": 1.3, "mean_gradient": 5.0,
        "clinical_notes": "WSS drop, concentric remodeling"
    },
    {
        "scenario": "diabetes", 
        "HR": 78, "SBP": 135, "DBP": 80, "PWV": 9.3, 
        "compliance_change": -0.40, "pulse_pressure": 55, "CO": 5.0, "EF": 55,
        "IMT": 0.85, "AVA": 3.5, "peak_velocity": 1.3, "mean_gradient": 5.0,
        "clinical_notes": "Glycocalyx damage, stiffness"
    },
    {
        "scenario": "aortic_stenosis", 
        "HR": 78, "SBP": 120, "DBP": 80, "PWV": 7.0, 
        "compliance_change": 0.0, "pulse_pressure": 40, "CO": 4.5, "EF": 50,
        "IMT": 0.65, "AVA": 1.0, "peak_velocity": 3.4, "mean_gradient": 25.0,
        "clinical_notes": "Valve calcification"
    },
    {
        "scenario": "combined", 
        "HR": 82, "SBP": 150, "DBP": 80, "PWV": 10.5, 
        "compliance_change": -0.45, "pulse_pressure": 70, "CO": 5.0, "EF": 45,
        "IMT": 0.90, "AVA": 3.5, "peak_velocity": 1.4, "mean_gradient": 6.0,
        "clinical_notes": "AIx=35%, High risk (OR=13.0)"
    },
    
    # ОСТРЫЕ СОСТОЯНИЯ (СМЕРТЕЛЬНЫЕ ИСХОДЫ ИЗ CSV)
    {
        "scenario": "heart_attack", 
        "HR": 95, "SBP": 85, "DBP": 60, "PWV": 5.5, 
        "compliance_change": -0.20, "pulse_pressure": 25, "CO": 2.5, "EF": 25,
        "IMT": 0.65, "AVA": 3.5, "peak_velocity": 0.7, "mean_gradient": 2.0,
        "clinical_notes": "Cardiogenic Shock, Systolic Failure"
    },
    {
        "scenario": "stroke", 
        "HR": 110, "SBP": 190, "DBP": 110, "PWV": 12.0, 
        "compliance_change": -0.50, "pulse_pressure": 80, "CO": 5.5, "EF": 50,
        "IMT": 0.80, "AVA": 3.5, "peak_velocity": 1.6, "mean_gradient": 7.0,
        "clinical_notes": "HTN Crisis, AFib"
    },
    {
        "scenario": "cardiac_arrest", 
        "HR": 0, "SBP": 0, "DBP": 0, "PWV": 0.0, 
        "compliance_change": 0.0, "pulse_pressure": 0, "CO": 0.0, "EF": 0,
        "IMT": 0.65, "AVA": 3.5, "peak_velocity": 0.0, "mean_gradient": 0.0,
        "clinical_notes": "V-Fib / Asystole, CO=0"
    },
]

# ==============================================================================
# 2. ВСТРОЕННЫЕ ДАННЫЕ (FALLBACK) - ЭПИДЕМИОЛОГИЯ
# ==============================================================================
BASELINE_DATA = [
    {
        "condition": "Aortic thrombus",
        "icd10_code": "I74.0",
        "location": "Aorta",
        "mortality_pct": 25.0,
        "sex_ratio_M_F": 2.5,
        "peak_age_years": 68,
        "annual_incidence_per_100k": 3.2,
        "notes": "Saddle thrombus; high mortality if mesenteric involvement.",
        "source_url": "https://academic.oup.com/eurheartj/article/45/36/3538/7741105",
    },
    {
        "condition": "Pulmonary embolism",
        "icd10_code": "I26.9",
        "location": "Lung",
        "mortality_pct": 5.0,
        "sex_ratio_M_F": 1.2,
        "peak_age_years": 65,
        "annual_incidence_per_100k": 60.0,
        "notes": "Massive PE mortality 30-50%. Overall ~60/100k.",
        "source_url": "https://academic.oup.com/eurheartj/article/41/4/543/5556136",
    },
    {
        "condition": "Deep vein thrombosis",
        "icd10_code": "I80.2",
        "location": "Vein (leg)",
        "mortality_pct": 0.5,
        "sex_ratio_M_F": 1.1,
        "peak_age_years": 60,
        "annual_incidence_per_100k": 120.0,
        "notes": "3-month mortality mainly from PE.",
        "source_url": "https://www.ahajournals.org/doi/10.1161/CIR.0000000000000919",
    },
    {
        "condition": "Ischemic stroke (cardioembolic)",
        "icd10_code": "I63.4",
        "location": "Brain",
        "mortality_pct": 20.0,
        "sex_ratio_M_F": 1.1,
        "peak_age_years": 72,
        "annual_incidence_per_100k": 30.0,
        "notes": "Cardioembolic subtype highest mortality.",
        "source_url": "https://academic.oup.com/eurheartj/article/42/5/373/5899003",
    },
    {
        "condition": "Coronary thrombosis (STEMI)",
        "icd10_code": "I21.0",
        "location": "Heart",
        "mortality_pct": 7.0,
        "sex_ratio_M_F": 2.0,
        "peak_age_years": 62,
        "annual_incidence_per_100k": 80.0,
        "notes": "In-hospital mortality ~5-10% with PCI.",
        "source_url": "https://academic.oup.com/eurheartj/article/44/37/3720/7243210",
    },
    {
        "condition": "Mesenteric artery thrombosis",
        "icd10_code": "K55.0",
        "location": "Mesentery",
        "mortality_pct": 50.0,
        "sex_ratio_M_F": 1.3,
        "peak_age_years": 70,
        "annual_incidence_per_100k": 1.5,
        "notes": "Catastrophic; often diagnosed late.",
        "source_url": "https://www.wjgnet.com/1007-9327/full/v26/i29/3928.htm",
    },
    {
        "condition": "Renal artery thrombosis",
        "icd10_code": "I28.0",
        "location": "Kidney",
        "mortality_pct": 10.0,
        "sex_ratio_M_F": 1.4,
        "peak_age_years": 65,
        "annual_incidence_per_100k": 2.0,
        "notes": "Can lead to renal failure.",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/?term=renal+artery+thrombosis+incidence",
    },
    {
        "condition": "Portal vein thrombosis",
        "icd10_code": "I81",
        "location": "Liver/Portal",
        "mortality_pct": 5.0,
        "sex_ratio_M_F": 1.6,
        "peak_age_years": 55,
        "annual_incidence_per_100k": 4.0,
        "notes": "Cirrhosis-related in 25%.",
        "source_url": "https://www.journal-of-hepatology.eu/article/S0168-8278(22)00015-0/fulltext",
    },
    {
        "condition": "Carotid artery thrombosis",
        "icd10_code": "I63.2",
        "location": "Carotid",
        "mortality_pct": 15.0,
        "sex_ratio_M_F": 1.8,
        "peak_age_years": 68,
        "annual_incidence_per_100k": 5.0,
        "notes": "Often embolic origin; TIA precursor.",
        "source_url": "https://www.escardio.org/Guidelines/Clinical-Practice-Guidelines/2021-ESC-Guidelines-on-cardiovascular-disease-prevention-in-clinical-practice",
    },
    {
        "condition": "Cerebral venous sinus thrombosis",
        "icd10_code": "G08",
        "location": "Brain (sinus)",
        "mortality_pct": 5.0,
        "sex_ratio_M_F": 0.4,
        "peak_age_years": 40,
        "annual_incidence_per_100k": 1.5,
        "notes": "More common in young women (OCP).",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/28393945/",
    },
    {
        "condition": "Subclavian vein thrombosis",
        "icd10_code": "I82.8",
        "location": "Upper limb vein",
        "mortality_pct": 1.0,
        "sex_ratio_M_F": 1.3,
        "peak_age_years": 45,
        "annual_incidence_per_100k": 2.0,
        "notes": "Paget-Schroetter syndrome; effort thrombosis.",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/22459779/",
    },
    {
        "condition": "Atrial thrombus (LAA)",
        "icd10_code": "I51.3",
        "location": "Heart (LA)",
        "mortality_pct": 8.0,
        "sex_ratio_M_F": 1.2,
        "peak_age_years": 70,
        "annual_incidence_per_100k": 10.0,
        "notes": "Main source of cardioembolic stroke in AF.",
        "source_url": "https://www.escardio.org/Guidelines/Clinical-Practice-Guidelines/Atrial-Fibrillation-Management",
    },
]

# ==============================================================================
# 3. БИОФИЗИЧЕСКИЕ ПАРАМЕТРЫ ДЛЯ CFD (ТРИГГЕРЫ ТРОМБОЗА)
# ==============================================================================
BIOPHYSICAL_DATA = [
    {
        "trigger_name": "Low WSS",
        "clinical_condition": "Atherosclerosis",
        "cfd_critical_value": "WSS < 0.4 - 0.5 Pa",
        "mechanism": "Локальная стагнация потока. ЛПНП проникает в стенку.",
        "cfd_implementation": "Уравнение адвекции-диффузии. Оседание при WSS < 0.4 Па.",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/10591296/",
        "paper_title": "Hemodynamic shear stress and its role in atherosclerosis (JAMA, 1999)"
    },
    {
        "trigger_name": "High WSS",
        "clinical_condition": "Plaque Rupture & SIPA",
        "cfd_critical_value": "WSS > 50.0 - 100.0 Pa",
        "mechanism": "Разворачивает фактор фон Виллебранда (vWF). Активация тромбоцитов.",
        "cfd_implementation": "Lagrangian Particle Tracking. WSS > 50 Па активирует частицу.",
        "source_url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4312743/",
        "paper_title": "Shear stress-induced platelet activation (2015)"
    },
    {
        "trigger_name": "Hyperglycemia",
        "clinical_condition": "Endothelial dysfunction",
        "cfd_critical_value": "Glucose > 7.0 mmol/L",
        "mechanism": "Высокий сахар разрушает гликокаликс. Стенка становится сверх-липкой.",
        "cfd_implementation": "Снизить Compliance. Увеличить коэффициент прилипания.",
        "source_url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2805904/",
        "paper_title": "Diabetes and Atherothrombosis (Circulation, 2010)"
    },
    {
        "trigger_name": "Hypercholesterolemia",
        "clinical_condition": "Dyslipidaemia",
        "cfd_critical_value": "LDL > 3.0 mmol/L",
        "mechanism": "Физическое изменение геометрии (стеноз).",
        "cfd_implementation": "Mesh deformation — сужение просвета в областях с Low WSS.",
        "source_url": "https://academic.oup.com/eurheartj/article/41/1/111/5556353",
        "paper_title": "ESC/EAS Guidelines dyslipidaemias (2019)"
    },
    {
        "trigger_name": "Stent Struts",
        "clinical_condition": "In-stent Thrombosis",
        "cfd_critical_value": "Strut thickness 80 - 150 um",
        "mechanism": "Балки стента создают зоны отрыва потока и рециркуляции.",
        "cfd_implementation": "Добавление микро-препятствий (балок) в CFD-домен.",
        "source_url": "https://pubmed.ncbi.nlm.nih.gov/24585746/",
        "paper_title": "Impact of stent strut thickness on thrombosis (2014)"
    }
]

# ==============================================================================
# 4. ГЛОБАЛЬНЫЕ ПРИЧИНЫ РАЗРУШЕНИЯ АРТЕРИЙ (WHO / GBD DATA)
# ==============================================================================
WHO_CVD_RISK_FACTORS = [
    {
        "risk_factor_ru": "Высокое давление (Гипертония)",
        "mechanics": "Механическое перенапряжение стенки (High Wall Tension)",
        "attributable_cvd_deaths_millions": 10.8,
        "url": "https://www.who.int/news-room/fact-sheets/detail/hypertension"
    },
    {
        "risk_factor_ru": "Высокий холестерин (ЛПНП)",
        "mechanics": "Формирование атеросклеротической бляшки (Low WSS зоны)",
        "attributable_cvd_deaths_millions": 4.4,
        "url": "https://www.who.int/data/gho/indicator-metadata-registry/imr-details/3236"
    },
    {
        "risk_factor_ru": "Высокий сахар (Диабет)",
        "mechanics": "Разрушение гликокаликса, кальцификация, потеря эластичности",
        "attributable_cvd_deaths_millions": 3.2,
        "url": "https://www.who.int/news-room/fact-sheets/detail/diabetes"
    },
    {
        "risk_factor_ru": "Курение",
        "mechanics": "Эндотелиальная дисфункция, системное воспаление",
        "attributable_cvd_deaths_millions": 3.3,
        "url": "https://www.who.int/news-room/fact-sheets/detail/tobacco"
    },
    {
        "risk_factor_ru": "Ожирение (Высокий ИМТ)",
        "mechanics": "Системное воспаление, рост сердечного выброса",
        "attributable_cvd_deaths_millions": 2.5,
        "url": "https://www.who.int/news-room/fact-sheets/detail/obesity-and-overweight"
    },
]

# ==============================================================================
# 5. WHO DATA API / KAGGLE (оставлено для целостности изначального скрипта)
# ==============================================================================
def fetch_who_data() -> pd.DataFrame | None:
    API_URL  = "https://ghoapi.azureedge.net/api/NCDMORT3070"
    try:
        r = requests.get(API_URL, params={"$top": "2000", "$filter": "Dim1 eq 'BTSX'"}, timeout=20)
        r.raise_for_status()
        records = r.json().get("value", [])
        if not records: return None
        df = pd.DataFrame(records)
        keep = [c for c in ["SpatialDim", "SpatialDimType", "TimeDim", "Dim1", "NumericValue"] if c in df.columns]
        df = df[keep].rename(columns={"SpatialDim": "country_code", "TimeDim": "year", "NumericValue": "prob_death_30_70_pct"})
        df.to_csv(WHO_OUT, index=False, encoding="utf-8-sig")
        return df
    except: return None

# ==============================================================================
# 6. ГЕНЕРАЦИЯ PDF С ГРАФИКАМИ, ТАБЛИЦЕЙ И ССЫЛКАМИ
# ==============================================================================
def generate_pdf_report(df: pd.DataFrame, risk_df: pd.DataFrame, hemo_df: pd.DataFrame, bio_df: pd.DataFrame, who_df: pd.DataFrame | None):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return

    pdf_path = OUT_DIR / "thrombosis_report.pdf"
    print(f"\n[PDF] Генерация отчета: {pdf_path}...")

    with PdfPages(pdf_path) as pdf:
        # Страница 1: Летальность vs Встречаемость
        fig, ax = plt.subplots(figsize=(12, 7))
        colors = ['red' if m > 15 else 'orange' if m > 5 else 'blue' for m in df['mortality_pct']]
        sizes = df['peak_age_years'] * 5
        ax.scatter(df['annual_incidence_per_100k'], df['mortality_pct'], s=sizes, c=colors, alpha=0.6, edgecolors='k')
        for i, row in df.iterrows():
            ax.annotate(row['condition'], (row['annual_incidence_per_100k'], row['mortality_pct']), xytext=(5, 5), textcoords='offset points', fontsize=8)
        ax.set_title("Летальность vs Встречаемость тромбозов (размер = возраст группы)", fontsize=14)
        ax.set_xlabel("Случаев на 100 000 человек в год"); ax.set_ylabel("Летальность (%)")
        ax.grid(True, linestyle='--', alpha=0.7)
        pdf.savefig(fig, bbox_inches='tight'); plt.close()

        # Страница 2: От чего разрушаются артерии
        fig, ax = plt.subplots(figsize=(12, 7))
        risk_sorted = risk_df.sort_values("attributable_cvd_deaths_millions", ascending=True)
        bars = ax.barh(range(len(risk_sorted)), risk_sorted['attributable_cvd_deaths_millions'], color='darkred', alpha=0.8)
        ax.set_yticks(range(len(risk_sorted))); ax.set_yticklabels(risk_sorted['risk_factor_ru'], fontsize=11)
        ax.set_title("От чего разрушается артерия: Главные факторы риска (ВОЗ / GBD)", fontsize=14)
        ax.set_xlabel("Смертей от ССЗ в мире ежегодно (миллионы человек)", fontsize=12)
        ax.grid(axis='x', linestyle='--', alpha=0.7)
        for bar in bars:
            width = bar.get_width()
            ax.annotate(f'{width} млн', xy=(width, bar.get_y() + bar.get_height() / 2), xytext=(3, 0), textcoords="offset points", ha='left', va='center', fontsize=11, fontweight='bold')
        mech_text = "\n".join([f"• {r['risk_factor_ru']}:\n  {r['mechanics']}" for _, r in risk_df.iterrows()])
        plt.figtext(0.1, -0.2, f"Как именно фактор портит артерию:\n{mech_text}", fontsize=10, ha="left", bbox={"facecolor":"lightgrey", "alpha":0.5, "pad":5})
        pdf.savefig(fig, bbox_inches='tight'); plt.close()

        # Страница 3: Целевые Клинические Показатели
        fig, ax = plt.subplots(figsize=(16, 7))
        # Добавим внутренние отступы (margin), чтобы таблица не упиралась в правый край PDF
        fig.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.05)
        
        ax.axis('off')
        
        # Форматируем данные для таблицы
        table_data = []
        for _, r in hemo_df.iterrows():
            sbp = r['SBP'] if pd.notnull(r['SBP']) else '-'
            dbp = r['DBP'] if pd.notnull(r['DBP']) else '-'
            pwv = r['PWV'] if pd.notnull(r['PWV']) else '-'
            comp = f"{r['compliance_change']*100:.0f}%" if pd.notnull(r['compliance_change']) else '-'
            
            # Собираем маркеры в строку
            markers = []
            if pd.notnull(r['CO']): markers.append(f"CO={r['CO']}")
            if pd.notnull(r['EF']): markers.append(f"EF={r['EF']}%")
            if pd.notnull(r['IMT']): markers.append(f"IMT={r['IMT']}")
            if pd.notnull(r['AVA']): markers.append(f"AVA={r['AVA']}")
            if pd.notnull(r['mean_gradient']): markers.append(f"Grad={r['mean_gradient']}")
            
            notes = r.get('clinical_notes', '')
            if notes:
                markers.append(notes)
            
            table_data.append([
                r['scenario'], 
                r['HR'], 
                f"{sbp} / {dbp}", 
                pwv, 
                comp, 
                ", ".join(markers)
            ])
            
        # Делаем последнюю колонку еще шире
        col_widths = [0.14, 0.06, 0.10, 0.10, 0.12, 0.48]
        
        tbl = ax.table(cellText=table_data, colLabels=["Сценарий", "HR", "SBP / DBP", "PWV (m/s)", "Δ Compliance", "Специфические маркеры"], 
                       cellLoc='left', loc='center', colWidths=col_widths)
        tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1.0, 2.5)
        
        for (r, c), cell in tbl.get_celld().items():
            cell.PAD = 0.03
            if r == 0: 
                cell.set_facecolor('#4472C4')
                cell.set_text_props(color='white', fontweight='bold')
            elif r % 2 == 0: 
                cell.set_facecolor('#f0f0f0')
                
        fig.suptitle("Клинические таргеты для 0D-двигателя", fontsize=16, fontweight='bold')
        pdf.savefig(fig, bbox_inches='tight'); plt.close()

        # Страница 4: Источники и ссылки
        fig, ax = plt.subplots(figsize=(12, 11))
        ax.axis('off')
        fig.suptitle("Источники данных и ссылки на исследования", fontsize=16, fontweight='bold', y=0.95)
        
        y_pos = 0.90
        def add_text(text, fontsize, bold=False, color='black'):
            nonlocal y_pos
            fw = 'bold' if bold else 'normal'
            fig.text(0.05, y_pos, text, fontsize=fontsize, fontweight=fw, color=color, wrap=True)
            y_pos -= 0.02

        add_text("ЭПИДЕМИОЛОГИЯ И ЛЕТАЛЬНОСТЬ (Clinical Guidelines):", 12, True, 'darkblue'); y_pos -= 0.01
        for _, row in df.iterrows():
            add_text(f"• {row['condition']}: {row['source_url']}", 9)
            if y_pos < 0.05: break
        
        y_pos -= 0.03
        add_text("БИОФИЗИЧЕСКИЕ ТРИГГЕРЫ ТРОМБОЗОВ (CFD Mechanics):", 12, True, 'darkblue'); y_pos -= 0.01
        for _, row in bio_df.iterrows():
            add_text(f"• {row['clinical_condition']} ({row['paper_title']}): {row['source_url']}", 9)

        y_pos -= 0.03
        add_text("ГЛОБАЛЬНЫЕ ФАКТОРЫ РИСКА (WHO):", 12, True, 'darkblue'); y_pos -= 0.01
        for _, row in risk_df.iterrows():
            add_text(f"• {row['risk_factor_ru']}: {row['url']}", 9)
            
        pdf.savefig(fig, bbox_inches='tight'); plt.close()

# ==============================================================================
# ТОЧКА ВХОДА
# ==============================================================================
def main():
    print("=" * 70)
    print("THROMBOSIS & CLINICAL TARGETS DATA FETCHER")
    print(f"Вывод в: {OUT_DIR}")
    print("=" * 70)

    df = pd.DataFrame(BASELINE_DATA)
    # Маппинг эпидемиологии на 0D-сценарии
    df["model_scenario"] = df["condition"].map(CAUSE_TO_SCENARIO).fillna("baseline")

    bio_df = pd.DataFrame(BIOPHYSICAL_DATA)
    risk_df = pd.DataFrame(WHO_CVD_RISK_FACTORS)
    hemo_df = pd.DataFrame(CLINICAL_TARGETS)

    who_df = None
    if SOURCE in ("all", "who"):
        who_df = fetch_who_data()

    df.to_csv(CSV_OUT, index=False, encoding="utf-8-sig")
    bio_df.to_csv(BIOPHY_OUT, index=False, encoding="utf-8-sig")
    risk_df.to_csv(RISK_OUT, index=False, encoding="utf-8-sig")
    hemo_df.to_csv(HEMO_OUT, index=False, encoding="utf-8-sig")

    generate_pdf_report(df, risk_df, hemo_df, bio_df, who_df)
    
    print("\n[DONE] Готовые файлы:")
    print(f"  {CSV_OUT}")
    print(f"  {HEMO_OUT}  <-- Гемодинамика для 0D")
    print(f"  {BIOPHY_OUT}")
    print(f"  {RISK_OUT}")
    print(f"  {OUT_DIR / 'thrombosis_report.pdf'}")

if __name__ == "__main__":
    main()
