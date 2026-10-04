# -*- coding: utf-8 -*-
"""
================================================================================
HEART ENGINE SWEEP v4 — 0D гемодинамическая модель, 5 клинических сценариев.
================================================================================
Калибровка: normal(120/80) → mild → moderate → severe → severe combined.
Исправлен расчёт CO, стартовые объёмы приведены в равновесие.
Включены точные цифры из клинических исследований.
================================================================================
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from dataclasses import dataclass
from copy import deepcopy
import sys, os

if sys.platform != "win32" and os.path.exists("/mnt/c/Work/Cursor_CFD"):
    BASE_DIR = "/mnt/c/Work/Cursor_CFD"
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))

OUT_DIR = os.path.join(BASE_DIR, "0D_results")
os.makedirs(OUT_DIR, exist_ok=True)

DT            = 0.0005
WARMUP_CYCLES = 25        # 25 циклов для идеального равновесия
RECORD_CYCLES = 5
D_INLET_MM    = 20.0

# ==============================================================================
# БАЗОВЫЕ ПАРАМЕТРЫ  —  здоровый человек, 70 кг, покой (120/80, CO~5 L/min)
# ==============================================================================
BASE_CHAMBERS = {
    "LV": dict(Emax=2.00, Emin=0.060, V0=5,   start=0.16, Tmax=0.30),
    "RV": dict(Emax=0.55, Emin=0.040, V0=5,   start=0.16, Tmax=0.30),
    "LA": dict(Emax=0.15, Emin=0.070, V0=5,   start=0.00, Tmax=0.09),
    "RA": dict(Emax=0.12, Emin=0.050, V0=5,   start=0.00, Tmax=0.09),
}

BASE_VESSELS = {
    "AO": dict(C=1.45,  V0=250.0),    # аорта + крупные артерии
    "PA": dict(C=4.30,  V0=60.0),     # лёгочная артерия
    "SV": dict(C=60.0,  V0=2700.0),   # системные вены
    "PV": dict(C=15.0,  V0=300.0),    # лёгочные вены
}

BASE_R = dict(
    mitral=0.005,   tricuspid=0.005,
    aortic=0.008,   pulmonic=0.005,
    systemic=1.35,  pulmonary=0.08,
    sv_ra=0.05,     pv_la=0.05,
)

# Идеально подобранные стартовые объемы под 120/80 (отсутствие скачков на старте)
INIT_VOLUMES = {
    'LA': 45.0,  'LV': 135.0,
    'RA': 45.0,  'RV': 135.0,
    'AO': 380.0, 'PA': 110.0,
    'SV': 3100.0,'PV': 380.0,
}

ECG_BASE = {
    'P': dict(t=0.02,  amp=0.15,  w=0.040),
    'Q': dict(t=0.150, amp=-0.12, w=0.010),
    'R': dict(t=0.165, amp=1.00,  w=0.012),
    'S': dict(t=0.180, amp=-0.18, w=0.012),
    'T': dict(t=0.380, amp=0.28,  w=0.070),
}

SOUNDS = dict(
    S1=dict(freq=45, dur=0.14, decay=28),
    S2=dict(freq=90, dur=0.11, decay=35),
    S3=dict(freq=35, dur=0.10, decay=40),
    S4=dict(freq=28, dur=0.10, decay=40),
)

# ==============================================================================
# ТОЧНЫЕ КЛИНИЧЕСКИЕ ДАННЫЕ ИЗ НАУЧНЫХ СТАТЕЙ
# ==============================================================================
CLINICAL_DATA = {
    "Гипертония": [
        ("Aortic Distensibility", "3.5 ± 0.7", "1.4 ± 0.3 (-60%)", "10⁻⁶ cm²/dyne", "Stefanadis, Circulation 1997"),
        ("Carotid-Femoral PWV", "6.5 - 7.2", "8.8 - 11.4", "m/s", "Ref Values Collab, EHJ 2010"),
        ("Mean Aortic Wall Thick.", "2.23 (max 3.4)", "2.45 (max 3.6)", "mm", "Malayeri, Am J Cardiol 2008"),
        ("Peak Systolic WSS", "1.2 - 1.8", "0.3 - 0.8", "Pa", "Callaghan, J Magn Reson Imaging 2011"),
        ("Peak Tensile Wall Stress", "100 - 150", "250 - 380 (+150%)", "kPa", "Beller, J Thorac Cardiovasc Surg 2004"),
    ],
    "Диабет 2 типа": [
        ("Carotid-Femoral PWV", "8.0 ± 1.6", "9.3 ± 2.0 (+16%)", "m/s", "Laugesen, Diabetes Care 2013"),
        ("Glycocalyx Thickness", "0.9 ± 0.1", "0.5 ± 0.1 (-44%)", "µm", "Nieuwdorp, Diabetes 2006"),
        ("Systemic Glycocalyx Vol", "1.5 ± 0.1", "0.8 ± 0.4 (-47%)", "Liters", "Nieuwdorp, Diabetes 2006"),
        ("Aortic Compliance Drop", "0%", "-39% to -43%", "%", "Lee, Diab Vasc Dis Res 2007"),
        ("Carotid IMT", "0.72 ± 0.11", "0.85 ± 0.14", "mm", "Brohall, Diabetic Medicine 2006"),
    ],
    "Холестерин / Аортальный стеноз": [
        ("Aortic Valve Area (AVA)", "3.0 - 4.0", "1.0 - 1.5 (mod)", "cm²", "Otto, JACC 2021 ACC/AHA"),
        ("Peak Aortic Velocity", "1.0 - 1.7", "3.0 - 3.9 (mod)", "m/s", "Otto, JACC 2021 ACC/AHA"),
        ("Mean Transvalvular Grad", "< 5.0", "20 - 39 (mod)", "mmHg", "Otto, JACC 2021 ACC/AHA"),
        ("LV Wall Thickness (LVH)", "6.0 - 10.0", "12.0 - 16.0+", "mm", "Dweck, JACC 2012"),
        ("Arch Atheroma Plaque", "< 2.0", "≥ 4.0 (OR 9.1 stroke)", "mm", "Amarenco, NEJM 1994"),
    ],
    "Комбинированная патология": [
        ("MI Odds Ratio (HTN)", "1.0", "1.91", "OR", "Yusuf, INTERHEART Lancet 2004"),
        ("MI Odds Ratio (HTN+DM+Sm)", "1.0", "13.01 (PAR 53%)", "OR", "Yusuf, INTERHEART Lancet 2004"),
        ("Annual PWV Increase", "5.7 ± 1.1", "68.3 ± 7.1", "cm/s/yr", "Tomiyama, Hypertension 2006"),
        ("Central Pulse Pressure", "30.0 - 40.0", "60.0 - 78.0+", "mmHg", "Prenner & Chirinos, Atherosclerosis 2015"),
        ("Augmentation Index", "5.0 - 15.0", "30.0 - 45.0+", "%", "Chirinos, ATVB 2019"),
    ],
}

# ==============================================================================
# КЛИНИЧЕСКИЕ СЦЕНАРИИ (Хронические + Острые исходы из CSV)
# ==============================================================================
SCENARIOS = [
    # 0. ЗДОРОВЫЙ
    dict(name="0_baseline",
         HR=72.0, R_aortic=0.008, Emax_LV=2.00,
         Emax_LA=0.15, Emin_LA=0.070,
         C_sys_factor=1.00, R_sys_factor=1.00,
         desc="Здоровый (Baseline)"),

    # 1. ГИПЕРТОНИЯ
    dict(name="1_hypertension_mild",
         HR=76.0, R_aortic=0.008, Emax_LV=2.10,
         Emax_LA=0.15, Emin_LA=0.070,
         C_sys_factor=0.90, R_sys_factor=1.25,
         desc="Гипертония (↑R_sys, ↓C)"),

    # 2. ДИАБЕТ / ЖЁСТКОСТЬ
    dict(name="2_diabetes_stiff",
         HR=78.0, R_aortic=0.008, Emax_LV=2.05,
         Emax_LA=0.15, Emin_LA=0.070,
         C_sys_factor=0.60, R_sys_factor=1.05,
         desc="Диабет: жёсткость (↓C 40%)"),

    # 3. АОРТАЛЬНЫЙ СТЕНОЗ
    dict(name="3_aortic_stenosis",
         HR=78.0, R_aortic=0.065, Emax_LV=2.30,
         Emax_LA=0.15, Emin_LA=0.070,
         C_sys_factor=1.00, R_sys_factor=1.00,
         desc="Аортальный стеноз (R_av ×8)"),

    # 4. ТЯЖЁЛАЯ КОМБИ (Только сосудистые патологии)
    dict(name="4_severe_combined",
         HR=82.0, R_aortic=0.008, Emax_LV=2.10,
         Emax_LA=0.15, Emin_LA=0.070,
         C_sys_factor=0.55, R_sys_factor=1.25,
         desc="Комбо: АГ + Диабет (Резко ↓C, ↑R_sys)"),

    # ================== ОСТРЫЕ ИСХОДЫ (МАППИНГ ИЗ CSV) ==================

    # 5. ИНФАРКТ МИОКАРДА (Heart Attack)
    dict(name="5_heart_attack",
         HR=95.0, R_aortic=0.008, Emax_LV=0.90,   # Резкое падение сократимости
         Emax_LA=0.15, Emin_LA=0.150,             # Диастолическая дисфункция ЛЖ/ЛП
         C_sys_factor=1.00, R_sys_factor=1.10,
         desc="Инфаркт миокарда (↓Emax, Кардиогенный шок)"),

    # 6. ИНСУЛЬТ (Stroke) - Гипертонический криз + Мерцательная аритмия
    dict(name="6_stroke",
         HR=110.0, R_aortic=0.008, Emax_LV=2.20,
         Emax_LA=0.05, Emin_LA=0.050,             # Мерцательная аритмия (нет систолы предсердий)
         C_sys_factor=0.50, R_sys_factor=1.80,    # Гипертонический криз
         desc="Инсульт (Гипертонический криз + Фибрилляция предсердий)"),

    # 7. ОСТАНОВКА СЕРДЦА (Cardiac Arrest / Sudden Death)
    dict(name="7_cardiac_arrest",
         HR=0.0, R_aortic=0.008, Emax_LV=0.06,    # Emax = Emin (абсолютно плоская эластанса, нет сокращений)
         Emax_LA=0.07, Emin_LA=0.070,             # Предсердия тоже не сокращаются
         C_sys_factor=1.00, R_sys_factor=1.00,
         desc="Остановка сердца (Asystole, HR=0, CO=0)"),
]

# ==============================================================================
# ФИЗИКА
# ==============================================================================
def _elastance(tc, start, Tmax, Emax, Emin, cycle):
    tn = tc - start
    if tn < 0: tn += cycle
    if tn < Tmax: e = 0.5 * (1 - np.cos(np.pi * tn / Tmax))
    elif tn < 1.5 * Tmax: e = 0.5 * (1 + np.cos(2 * np.pi * (tn - Tmax) / Tmax))
    else: e = 0.0
    return Emin + (Emax - Emin) * e

@dataclass
class ValveState:
    is_open: bool = False
    flow: float   = 0.0

def _update_valve(p_up, p_dn, resistance):
    if p_up > p_dn: return ValveState(True, (p_up - p_dn) / resistance)
    return ValveState()

def _make_ecg_waves(scenario):
    waves = deepcopy(ECG_BASE)
    ratio = scenario["Emax_LV"] / 2.0
    waves['R']['amp'] *= ratio ** 0.4
    waves['S']['amp'] *= ratio ** 0.3
    if ratio > 1.1: waves['T']['amp'] *= max(0.6, 1.0 - 0.3 * (ratio - 1.1))
    return waves

def _synthesize_pcg(t_array, events):
    pcg = np.zeros_like(t_array)
    rng = np.random.RandomState(42)
    for te, name, scale in events:
        if name in SOUNDS:
            p = SOUNDS[name]
            mask = (t_array >= te) & (t_array <= te + p['dur'])
            if not np.any(mask): continue
            tau = t_array[mask] - te
            pcg[mask] += scale * np.exp(-p['decay'] * tau) * np.sin(2*np.pi*p['freq']*tau)
        elif name == 'MURMUR_SYS':
            dur = 0.18
            mask = (t_array >= te) & (t_array <= te + dur)
            if not np.any(mask): continue
            tau = t_array[mask] - te
            env = np.sin(np.pi * tau / dur)
            pcg[mask] += scale * env * rng.randn(len(tau)) * 0.25
    return pcg

# ==============================================================================
# СЕРДЦЕ
# ==============================================================================
class Heart:
    def __init__(self, scenario):
        self.scenario_name = scenario["name"]
        self.HR    = scenario["HR"]
        self.CYCLE = 60.0 / self.HR if self.HR > 0 else 1.0
        self.ecg_waves = _make_ecg_waves(scenario)

        self.chambers = deepcopy(BASE_CHAMBERS)
        self.chambers["LV"]["Emax"] = scenario["Emax_LV"]
        self.chambers["LA"]["Emax"] = scenario["Emax_LA"]
        self.chambers["LA"]["Emin"] = scenario["Emin_LA"]

        self.vessels = deepcopy(BASE_VESSELS)
        self.vessels["AO"]["C"] *= scenario.get("C_sys_factor", 1.0)
        self.vessels["SV"]["C"] *= scenario.get("C_sys_factor", 1.0)

        self.R = deepcopy(BASE_R)
        self.R["aortic"]    = scenario["R_aortic"]
        self.R["systemic"] *= scenario.get("R_sys_factor", 1.0)

        self.V = deepcopy(INIT_VOLUMES)
        self.valves       = {k: False for k in ['mitral','tricuspid','aortic','pulmonic']}
        self._prev_valves = dict(self.valves)
        self.time          = 0.0
        self.phase         = 1
        self.ecg           = 0.0
        self.venous_pulse  = 0.0
        self._prev_RV_P    = 0.0
        self.events        = []
        
        self._prev_cycle_id = 0
        self._lv_max = self.V["LV"]
        self._lv_min = self.V["LV"]
        self.last_EDV       = self.V["LV"]
        self.last_ESV       = self.V["LV"]
        self.stroke_volume  = 0.0
        self.cardiac_output = 0.0
        self.Q_AO           = 0.0
        self.P              = {}

        self._max_grad_this  = 0.0
        self._max_grad_prev  = 0.0

    def _pressures(self, tc):
        P = {}
        for n, c in self.chambers.items():
            P[n] = _elastance(tc, c['start'], c['Tmax'], c['Emax'], c['Emin'], self.CYCLE) * (self.V[n] - c['V0'])
        for n, v in self.vessels.items():
            P[n] = (self.V[n] - v['V0']) / v['C']
        return P

    def _update_phase(self, tc, P):
        ao, mi = self.valves['aortic'], self.valves['mitral']
        if tc < 1.5 * self.chambers['LA']['Tmax'] and mi and not ao:
            self.phase = 1; return
        if not mi and not ao:
            self.phase = 2 if self.V['LV'] >= self._lv_max - 1e-6 else 5; return
        if ao:
            f = (self.V['LV'] - self.last_ESV) / max(self.last_EDV - self.last_ESV, 1e-6)
            self.phase = 3 if f > 0.5 else 4; return
        if mi:
            f = (self.V['LV'] - self.last_ESV) / max(self.last_EDV - self.last_ESV, 1e-6)
            self.phase = 6 if f < 0.75 else 7; return
        self.phase = 7

    def step(self, dt=DT):
        # ---------------------------------------------------------------
        # ОСТАНОВКА СЕРДЦА: только пассивное затухание Windkessel.
        # Никаких клапанов, фаз, циклов — давление экспоненциально → 0.
        # ---------------------------------------------------------------
        if self.HR == 0:
            P = {}
            for n, c in self.chambers.items():
                P[n] = c['Emin'] * (self.V[n] - c['V0'])
            for n, v in self.vessels.items():
                P[n] = (self.V[n] - v['V0']) / v['C']

            # Только пассивные потоки через сосудистое сопротивление (нет насоса)
            q_sys = max(P["AO"] - P["SV"], 0) / self.R["systemic"]
            q_pul = max(P["PA"] - P["PV"], 0) / self.R["pulmonary"]
            q_svr = max(P["SV"] - P["RA"], 0) / self.R["sv_ra"]
            q_pvl = max(P["PV"] - P["LA"], 0) / self.R["pv_la"]

            self.V["AO"] -= dt * q_sys
            self.V["SV"] += dt * (q_sys - q_svr)
            self.V["PA"] -= dt * q_pul
            self.V["PV"] += dt * (q_pul - q_pvl)

            self.Q_AO = 0.0
            self.ecg  = 0.02 * np.sin(2 * np.pi * 45 * self.time) * np.exp(-self.time * 0.3)
            self.P    = P
            self.cardiac_output = 0.0
            self.stroke_volume  = 0.0
            self.time += dt
            return P

        tc = self.time % self.CYCLE
        cycle_id = int(self.time // self.CYCLE)

        if cycle_id != self._prev_cycle_id:
            self.last_EDV       = self._lv_max
            self.last_ESV       = self._lv_min
            self.stroke_volume  = self.last_EDV - self.last_ESV
            self.cardiac_output = self.stroke_volume * self.HR / 1000.0 # л/мин
            self._lv_max = self.V["LV"]
            self._lv_min = self.V["LV"]
            self._max_grad_prev = self._max_grad_this
            self._max_grad_this = 0.0
            self._prev_cycle_id = cycle_id

        self._lv_max = max(self._lv_max, self.V["LV"])
        self._lv_min = min(self._lv_min, self.V["LV"])

        if self.HR == 0:
            # Асистолия: нет пиков, только артефактный шум
            self.ecg = 0.02 * np.sin(2 * np.pi * 45 * self.time) * np.exp(-0.1 * self.time % 1.0)
        else:
            v_ecg = 0.0
            for w in self.ecg_waves.values():
                d = tc - w['t']
                v_ecg += w['amp'] * np.exp(-d*d / (2.0 * w['w']**2))
            self.ecg = v_ecg

        P = self._pressures(tc)

        mi = _update_valve(P["LA"], P["LV"], self.R["mitral"])
        ao = _update_valve(P["LV"], P["AO"], self.R["aortic"])
        tr = _update_valve(P["RA"], P["RV"], self.R["tricuspid"])
        pu = _update_valve(P["RV"], P["PA"], self.R["pulmonic"])
        self.valves = {'mitral': mi.is_open, 'aortic': ao.is_open, 'tricuspid': tr.is_open, 'pulmonic': pu.is_open}

        if ao.is_open:
            self._max_grad_this = max(self._max_grad_this, P["LV"] - P["AO"])

        pv = self._prev_valves
        if (pv['mitral'] and not mi.is_open) or (pv['tricuspid'] and not tr.is_open):
            if not any(e[1] == 'S1' and abs(e[0]-self.time) < 0.02 for e in self.events):
                self.events.append((self.time, 'S1', 1.0))
                murmur = min(self._max_grad_prev / 80.0, 1.0)
                if murmur > 0.08:
                    self.events.append((self.time + 0.04, 'MURMUR_SYS', murmur))
        if (pv['aortic'] and not ao.is_open) or (pv['pulmonic'] and not pu.is_open):
            if not any(e[1] == 'S2' and abs(e[0]-self.time) < 0.02 for e in self.events):
                s2 = max(0.3, 0.8 * (0.008 / self.R["aortic"]) ** 0.3)
                self.events.append((self.time, 'S2', s2))
        if (not pv['mitral']) and mi.is_open:
            self.events.append((self.time + 0.10, 'S3', 0.18))
        if abs(tc - self.chambers['LA']['Tmax']) < dt:
            self.events.append((self.time, 'S4', 0.10))
        self._prev_valves = dict(self.valves)

        q_sys = max(P["AO"]-P["SV"], 0) / self.R["systemic"]
        q_pul = max(P["PA"]-P["PV"], 0) / self.R["pulmonary"]
        q_svr = max(P["SV"]-P["RA"], 0) / self.R["sv_ra"]
        q_pvl = max(P["PV"]-P["LA"], 0) / self.R["pv_la"]

        self.V["LA"] += dt * (q_pvl - mi.flow)
        self.V["LV"] += dt * (mi.flow - ao.flow)
        self.V["RA"] += dt * (q_svr - tr.flow)
        self.V["RV"] += dt * (tr.flow - pu.flow)
        self.V["AO"] += dt * (ao.flow - q_sys)
        self.V["PA"] += dt * (pu.flow - q_pul)
        self.V["SV"] += dt * (q_sys - q_svr)
        self.V["PV"] += dt * (q_pul - q_pvl)

        dRV = (P['RV'] - self._prev_RV_P) / dt
        self.venous_pulse = P['RA'] + (0.010 * dRV if (not tr.is_open and dRV > 0) else 0.0)
        self._prev_RV_P = P['RV']

        self._update_phase(tc, P)
        self.Q_AO = ao.flow
        self.P    = P
        self.time += dt
        return P

# ==============================================================================
# ЗАПУСК
# ==============================================================================
def run_scenario(scenario, warmup=WARMUP_CYCLES, record=RECORD_CYCLES):
    heart = Heart(scenario)
    for _ in range(int(warmup * heart.CYCLE / DT)): heart.step()
    heart.events = []

    n = int(record * heart.CYCLE / DT)
    t0 = heart.time
    keys = ['ecg','LV','LA','RV','RA','AO','PA','V_LV','V_LA','V_RV','V_RA','phase','venous_pulse','Q_AO']
    data = {k: [] for k in keys}

    for _ in range(n):
        P = heart.step()
        data['ecg'].append(heart.ecg); data['phase'].append(heart.phase); data['venous_pulse'].append(heart.venous_pulse)
        for ch in ('LV','LA','RV','RA','AO','PA'): data[ch].append(P[ch])
        for ch in ('LV','LA','RV','RA'): data[f'V_{ch}'].append(heart.V[ch])
        data['Q_AO'].append(heart.Q_AO)

    for k in data: data[k] = np.array(data[k])
    t = np.arange(n) * DT + t0
    evs = [(te, nm, a) for te, nm, a in heart.events if t0 <= te <= t[-1] + 0.2]
    data['pcg'] = _synthesize_pcg(t, evs)

    # --- Статистика из последнего цикла ---
    CYCLE = heart.CYCLE
    idx_last = t >= (t[-1] - CYCLE)
    Q_last = data['Q_AO'][idx_last]
    P_ao_last = data['AO'][idx_last]
    V_lv_last = data['V_LV'][idx_last]
    P_lv_last = data['LV'][idx_last]

    EDV = float(V_lv_last.max()); ESV = float(V_lv_last.min()); SV = EDV - ESV
    EF  = 100.0 * SV / EDV if EDV > 0 else 0.0
    CO  = SV * heart.HR / 1000.0  # л/мин

    SBP = float(P_ao_last.max()); DBP = float(P_ao_last.min()); MAP = float(P_ao_last.mean())
    grad_inst = P_lv_last - P_ao_last
    peak_gradient = float(grad_inst.max())
    mean_gradient = float(np.mean(grad_inst[grad_inst > 0])) if np.any(grad_inst > 0) else 0.0

    return dict(
        scenario=scenario['name'], desc=scenario['desc'], HR=heart.HR, CYCLE=CYCLE,
        CO=CO, EF=EF, SV=SV, EDV=EDV, ESV=ESV, SBP=SBP, DBP=DBP, MAP=MAP,
        peak_gradient=peak_gradient, mean_gradient=mean_gradient, P_LV_max=float(P_lv_last.max()),
        t=t, data=data, Q_AO=data['Q_AO'], P_AO=data['AO'], V_LV=data['V_LV'], P_LV=data['LV'],
    )

def _plot_dashboard(ax_list, t, data, title):
    axes = ax_list
    axes[0].plot(t, data['ecg'], 'k', lw=0.8)
    axes[0].set_ylabel('ЭКГ (мВ)'); axes[0].set_title(title, fontsize=12, fontweight='bold')
    axes[1].plot(t, data['LV'], label='P_LV', color='crimson', lw=1.2)
    axes[1].plot(t, data['LA'], label='P_LA', color='orange', lw=0.8)
    axes[1].plot(t, data['AO'], label='P_AO', color='darkred', lw=1.0, ls='--')
    axes[1].set_ylabel('Давление\n(ммHg, левые)'); axes[1].legend(fontsize=7)
    axes[2].plot(t, data['RV'], label='P_RV', color='steelblue', lw=1.2)
    axes[2].plot(t, data['RA'], label='P_RA', color='deepskyblue', lw=0.8)
    axes[2].plot(t, data['PA'], label='P_PA', color='navy', lw=1.0, ls='--')
    axes[2].set_ylabel('Давление\n(ммHg, правые)'); axes[2].legend(fontsize=7)
    axes[3].plot(t, data['V_LV'], label='V_LV', color='crimson', lw=1.2)
    axes[3].plot(t, data['V_RV'], label='V_RV', color='steelblue', lw=1.0)
    axes[3].plot(t, data['V_LA'], label='V_LA', color='orange', lw=0.7)
    axes[3].plot(t, data['V_RA'], label='V_RA', color='deepskyblue', lw=0.7)
    axes[3].set_ylabel('Объём (мл)'); axes[3].legend(fontsize=7)
    axes[4].plot(t, data['venous_pulse'], color='teal', lw=0.9); axes[4].set_ylabel('Венозный\nпульс')
    axes[5].plot(t, data['pcg'], 'k', lw=0.5); axes[5].set_ylabel('ФКГ')
    axes[6].step(t, data['phase'], where='post', color='purple', lw=0.8)
    axes[6].set_ylabel('Фаза'); axes[6].set_xlabel('Время (с)'); axes[6].set_yticks(range(1, 8))
    for a in axes: a.grid(alpha=0.25)

if __name__ == "__main__":
    print("=" * 72)
    print("HEART ENGINE SWEEP v4 — Калибровка: 120/80 Baseline")
    print("=" * 72)

    results = []
    for scen in SCENARIOS:
        print(f"\n  [{scen['name']:28s}] HR={scen['HR']:4.0f} ...", end="", flush=True)
        res = run_scenario(scen)
        results.append(res)
        idx = res['t'] >= (res['t'][-1] - res['CYCLE'])
        np.savez(f"{OUT_DIR}/turbine_bc_{scen['name']}.npz", t=res['t'][idx]-res['t'][idx][0], 
                 Q_AO_m3_s=res['Q_AO'][idx]*1e-6, u_mean_m_s=res['Q_AO'][idx]*1e-6/(np.pi*((D_INLET_MM/2000.0)**2)),
                 P_AO_Pa=res['P_AO'][idx]*133.322, P_AO_mmHg=res['P_AO'][idx], Q_AO_mL_s=res['Q_AO'][idx], T_PERIOD=res['CYCLE'])
        print(f"  OK  ({res['SBP']:.0f}/{res['DBP']:.0f}, CO={res['CO']:.1f})")

    pdf_path = f"{OUT_DIR}/heart_engine_report.pdf"
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#9467bd', '#d62728']
    
    with PdfPages(pdf_path) as pdf:
        # 1. Clinical Data
        fig_ref = plt.figure(figsize=(16, 20))
        fig_ref.suptitle("Гемодинамические и структурные изменения аорты (Клинические данные)", fontsize=16, fontweight='bold', y=0.98)
        y = 0.93
        for cond, rows in CLINICAL_DATA.items():
            fig_ref.text(0.05, y, cond, fontsize=13, fontweight='bold', color='darkred'); y -= 0.025
            fig_ref.text(0.05, y, f"  {'Параметр':<25s} {'Норма':<20s} {'Патология':<25s} {'Ед.':<10s} Источник", fontsize=9, fontfamily='monospace', color='gray'); y -= 0.018
            for param, norm, path, unit, src in rows:
                fig_ref.text(0.05, y, f"  {param:<25s} {norm:<20s} {path:<25s} {unit:<10s} {src}", fontsize=9, fontfamily='monospace'); y -= 0.018
            y -= 0.015
        
        y -= 0.01; fig_ref.text(0.05, y, "Параметры симуляции:", fontsize=13, fontweight='bold', color='darkblue'); y -= 0.025
        fig_ref.text(0.05, y, f"  {'Сценарий':<25s} {'HR':>4s} {'R_av':>6s} {'Emax':>5s} {'C_f':>5s} {'R_f':>5s}   Описание", fontsize=9, fontfamily='monospace', color='gray'); y -= 0.018
        for sc in SCENARIOS:
            fig_ref.text(0.05, y, f"  {sc['name']:<25s} {sc['HR']:4.0f} {sc['R_aortic']:6.3f} {sc['Emax_LV']:5.2f} {sc['C_sys_factor']:5.2f} {sc['R_sys_factor']:5.2f}   {sc['desc']}", fontsize=9, fontfamily='monospace'); y -= 0.018
        pdf.savefig(fig_ref); plt.close(fig_ref)

        # 2. Sweep
        fig2, axes2 = plt.subplots(2, 2, figsize=(16, 12))
        ax_q, ax_p, ax_u, ax_pv = axes2[0,0], axes2[0,1], axes2[1,0], axes2[1,1]
        R_m = (D_INLET_MM/1000.0)/2.0; r_pts = np.linspace(-R_m, R_m, 200); n_pow = 7.0
        
        # Расширенный список цветов для всех 8 сценариев
        plot_colors = ['blue', 'orange', 'purple', 'green', 'red', 'black', 'brown', 'cyan']
        
        for i, res in enumerate(results):
            c = plot_colors[i % len(plot_colors)]; lbl = res['scenario']; idx = res['t'] >= (res['t'][-1] - res['CYCLE']); t_cyc = res['t'][idx] - res['t'][idx][0]
            ax_q.plot(t_cyc, res['Q_AO'][idx], label=lbl, color=c, lw=2); ax_p.plot(t_cyc, res['P_AO'][idx], label=lbl, color=c, lw=2)
            ax_pv.plot(res['V_LV'][idx], res['P_LV'][idx], label=lbl, color=c, lw=2)
            u_max = (np.max(res['Q_AO'][idx])*1e-6 / (np.pi*R_m**2)) * (n_pow+1)*(2*n_pow+1)/(2*n_pow**2)
            ax_u.plot(r_pts*1000, u_max*np.clip(1-np.abs(r_pts)/R_m, 0, None)**(1.0/n_pow), label=f"{lbl} (u_max={u_max:.2f} м/с)", color=c, lw=2)
        ax_q.set_title("Q_AO(t), мл/с"); ax_q.legend(fontsize=7); ax_p.set_title("P_AO(t), ммHg"); ax_p.legend(fontsize=7)
        ax_u.set_title("u(r)"); ax_u.legend(fontsize=6); ax_pv.set_title("PV-петля ЛЖ"); ax_pv.legend(fontsize=7)
        for a in [ax_q, ax_p, ax_u, ax_pv]: a.grid(alpha=0.4)
        pdf.savefig(fig2); plt.close(fig2)

        # 3. Dashboards (для ВСЕХ сценариев)
        for idx_s in range(len(results)):
            res = results[idx_s]
            fig_d, axes_d = plt.subplots(7, 1, figsize=(14, 22), sharex=True, gridspec_kw={'height_ratios': [1,1.5,1.5,1.2,1,0.8,0.8]})
            _plot_dashboard(axes_d, res['t'], res['data'], f"Dashboard: {SCENARIOS[idx_s]['desc']}")
            axes_d[-1].set_xlim(res['t'][0], res['t'][-1])
            fig_d.tight_layout(); pdf.savefig(fig_d); plt.close(fig_d)

    print(f"\nPDF сохранён: {pdf_path}\n")
    print(f"{'Scenario':25s} | {'HR':>3s} | {'SBP/DBP':>9s} | {'MAP':>5s} | {'CO':>5s} | {'EF':>5s} | {'SV':>4s} | {'Peak ΔP':>7s} | {'Mean ΔP':>7s}")
    print("-" * 95)
    for r in results:
        print(f"{r['scenario']:25s} | {r['HR']:3.0f} | {r['SBP']:3.0f}/{r['DBP']:<3.0f}   | {r['MAP']:5.1f} | {r['CO']:5.1f} | {r['EF']:5.1f} | {r['SV']:4.0f} | {r['peak_gradient']:7.1f} | {r['mean_gradient']:7.1f}")
    print("=" * 95)
