# CFD & Healthcare Simulation Roadmap
## Siemens Portfolio Project

---

## Timeline Overview

```mermaid
gantt
    title Blood Pump CFD & Healthcare Simulation Project Roadmap
    dateFormat YYYY-MM-DD
    axisFormat %b %Y
    
    section Phase Completion Status
    Phase 1-3: Meshing & Geometry       :done, phase1, 2024-01-01, 2024-06-30
    Phase 4: Hemodynamics & Aorta       :done, phase4, 2024-07-01, 2024-12-31
    Phase 5: Pathophysiology & 0D       :done, phase5, 2025-01-01, 2025-06-30
    Phase 6: FSI (Fluid-Structure)      :active, phase6, 2025-07-01, 2026-06-30
    Phase 7-9: Advanced Models & ML     :crit, phase789, 2026-07-01, 2027-12-31
    
    section Detailed Subtasks
    1.1: Mesher Switch                  :done, sub1_1, 2024-01-01, 30d
    
    2.1: Blade Comparison (N1,N2,N3)    :done, sub2_1, 2024-01-15, 60d
    2.2: Mesh Refinement Testing        :done, sub2_2, 2024-03-15, 45d
    2.3: Hemodynamics Analysis          :done, sub2_3, 2024-04-30, 30d
    
    3.1: Kármán Vortex Implementation   :done, sub3_1, 2024-05-15, 45d
    3.2: Unsteady Navier-Stokes         :done, sub3_2, 2024-06-01, 30d
    
    4.1: Aorta Geometry & Setup         :done, sub4_1, 2024-07-01, 45d
    4.2: DolfinX Integration            :done, sub4_2, 2024-08-01, 45d
    4.3: 0D Heart Model Connection      :done, sub4_3, 2024-09-01, 45d
    4.4: Windkessel & Outlets           :done, sub4_4, 2024-09-15, 30d
    4.5: Clot Modeling                  :done, sub4_5, 2024-10-01, 30d
    4.6: IPCS Fixes (Critical)          :done, sub4_6, 2024-10-15, 60d
    4.7: Scheme Comparison              :done, sub4_7, 2024-12-01, 30d
    
    5.1: Thrombosis Database            :done, sub5_1, 2025-01-01, 30d
    5.2: WHO/GBD Data Parser            :done, sub5_2, 2025-01-15, 45d
    5.3: CVD Risk Data Collection       :done, sub5_3, 2025-02-01, 45d
    5.4: Visualization & Report Gen     :done, sub5_4, 2025-03-01, 30d
    
    6.1: Heart Engine Physics           :done, sub6_1, 2025-03-15, 45d
    6.2: 0D Model (Windkessel)          :done, sub6_2, 2025-04-01, 30d
    6.3: 5 WHO Scenarios                :done, sub6_3, 2025-04-30, 60d
    6.4: NPZ Export & Visualization     :done, sub6_4, 2025-06-01, 30d
    
    7.1: 0D→2D CFD BC Integration       :crit, sub7_1, 2025-06-15, 60d
    7.2: WSS Field Comparison           :crit, sub7_2, 2025-08-01, 60d
    7.3: Womersley Profile (Optional)   :sub7_3, 2025-09-15, 45d
    7.4: Pulsatility Index Metrics      :sub7_4, 2025-10-01, 30d
    
    8.1: FEBio Setup                    :crit, sub8_1, 2025-07-01, 60d
    8.2: MFEM/MFEMiFSI Setup            :crit, sub8_2, 2025-08-15, 60d
    8.3: FSI Mesh Creation (2D)         :done, sub8_3, 2025-09-01, 30d
    8.4: BDF2-ALE FSI Solver            :active, sub8_4, 2025-09-15, 90d
    8.5: Benchmark 1 (2D Channel)       :sub8_5, 2025-12-01, 60d
    8.6: Benchmark 2 (Carotid)          :sub8_6, 2026-02-01, 60d
    8.7: Benchmark 3 (Aortic FSI)       :sub8_7, 2026-04-01, 90d
    8.8: Turbine Integration            :sub8_8, 2026-07-01, 60d
    8.9: Pulsatile Aorta + Turbine      :sub8_9, 2026-09-01, 90d
    8.10: Advanced Wall Model           :crit, sub8_10, 2026-11-01, 120d
    
    9.1: Global Circulation Model       :crit, sub9_1, 2026-07-01, 90d
    9.2: Plaque Modeling                :crit, sub9_2, 2026-10-01, 90d
    9.3: Rupture Risk Simulation        :crit, sub9_3, 2027-01-01, 90d
    9.4: Stent Placement                :sub9_4, 2027-04-01, 90d
    9.5: Hemolysis Analysis             :sub9_5, 2027-07-01, 60d
    
    10.1: ML Dataset (Kaggle)           :done, sub10_1, 2025-01-01, 30d
    10.2: Health Prediction Model       :crit, sub10_2, 2027-01-01, 120d
    10.3: Diagnostic Assistant (CV)     :crit, sub10_3, 2027-05-01, 120d
```

---

## Phase Breakdown & Status

### ✅ **Phase 1–3: Meshing & Geometry** (COMPLETED)
| Subtask | Status | Timeline |
|---------|--------|----------|
| Mesher: gmsh_from_step.py | ✅ DONE | Q1 2024 |
| Blade comparison (N1, N2, N3) | ✅ DONE | Q1–Q2 2024 |
| Mesh refinement study | ✅ DONE | Q2 2024 |
| Kármán vortex street | ✅ DONE | Q2 2024 |
| Unsteady Navier–Stokes | ✅ DONE | Q2 2024 |

### ✅ **Phase 4: Hemodynamics & Aorta** (COMPLETED)
| Subtask | Status | Timeline |
|---------|--------|----------|
| Aorta geometry & DolfinX | ✅ DONE | Q3–Q4 2024 |
| 0D heart model + Windkessel | ✅ DONE | Q3–Q4 2024 |
| Clot modeling & Re=500 | ✅ DONE | Q4 2024 |
| 3 scheme comparison (IPCS, Newton-CN, Newton-BDF2) | ✅ DONE | Q4 2024 |
| **Critical IPCS fixes** | ✅ DONE | Q4 2024 |
| — Fix A3 (velocity leakage) | ✅ DONE | Q4 2024 |
| — apply_lifting correction | ✅ DONE | Q4 2024 |
| — Skew-symmetric convection | ✅ DONE | Q4 2024 |
| — Artificial viscosity | ✅ DONE | Q4 2024 |
| — Backflow stabilization | ✅ DONE | Q4 2024 |
| — Adaptive dt + ABORT | ✅ DONE | Q4 2024 |
| — Windkessel sync | ✅ DONE | Q4 2024 |
| — Mass error < 1.5% | ✅ DONE | Q4 2024 |
| — Detailed diagnostics | ✅ DONE | Q4 2024 |

### ✅ **Phase 5: Pathophysiology & 0D Models** (COMPLETED)
| Subtask | Status | Timeline |
|---------|--------|----------|
| Thrombosis database (12 entities) | ✅ DONE | Q1 2025 |
| WHO/GBD data parser | ✅ DONE | Q1–Q2 2025 |
| CFD biophysical triggers (WSS, cholesterol) | ✅ DONE | Q2 2025 |
| CVD risk factor analysis | ✅ DONE | Q2 2025 |
| Thrombosis report visualization | ✅ DONE | Q2 2025 |
| **Heart engine sweep** | ✅ DONE | Q2–Q3 2025 |
| — Time-varying elastance | ✅ DONE | Q2 2025 |
| — 0D Windkessel model | ✅ DONE | Q2 2025 |
| — Valve diode logic | ✅ DONE | Q2 2025 |
| — 5 WHO scenarios | ✅ DONE | Q3 2025 |
| — NPZ export for CFD | ✅ DONE | Q3 2025 |
| — PV loops & velocity profiles | ✅ DONE | Q3 2025 |
| Power law inlet profile (1/7) | ✅ DONE | Q3 2025 |

### 🔄 **Phase 6: FSI (Fluid–Structure Interaction)** (IN PROGRESS)
| Subtask | Status | Timeline | Notes |
|---------|--------|----------|-------|
| Smoke-test budget (DT=2.5e-5 s, 400 steps, 15 min) | ✅ DONE | Q3 2025 | DolfinX/FEniCSx prototype validated |
| 2D conforming FSI mesh | ✅ DONE | Q3 2025 | `schematic_rects_fsi.msh` (fluid=100, wall=200, interface=6) |
| BDF2-ALE FSI solver | 🔄 IN PROGRESS | Q3–Q4 2025 | Traction transfer, elastic wall, mesh velocity, GCL |
| **Phase 7: 0D → 2D CFD Sync** | 🔴 NOT STARTED | Q2–Q3 2025 | BLOCKING: required before FSI coupling |
| — BC integration (Q_AO, P_AO) | 🔴 NOT DONE | Q2 2025 | Critical path item |
| — WSS field comparison (5 scenarios) | 🔴 NOT DONE | Q3 2025 | Depends on BC integration |
| — Womersley profile (optional) | ⏸️ OPTIONAL | Q3 2025 | Nice-to-have analytic alternative |
| — Pulsatility index (PI) | 🔴 NOT DONE | Q3 2025 | Hemodynamics metric |
| **Phase 8: Immersed FSI (MFEM + FEBio)** | 🔴 NOT STARTED | Q3 2025–Q2 2026 | **CRITICAL PATH** |
| — FEBio hyperelastic tissue setup | 🔴 NOT DONE | Q3 2025 | Dependency for FSI |
| — MFEM / MFEMiFSI integration | 🔴 NOT DONE | Q3–Q4 2025 | Production framework |
| — Benchmark 1: 2D pulsating channel | 🔴 NOT DONE | Q4 2025 | Validation milestone |
| — Benchmark 2: 3D carotid + deformable walls | 🔴 NOT DONE | Q1 2026 | Physiological geometry |
| — Benchmark 3: Aortic arch FSI | 🔴 NOT DONE | Q2 2026 | Full system validation |
| — Turbine + aorta integration | 🔴 NOT DONE | Q3 2026 | Blood pump in circuit |
| — Advanced wall model (nonlinear, spatially varying h) | 🔴 NOT DONE | Q3–Q4 2026 | Large deformation mechanics |
| — Balloon smoke test (arch expansion) | 🔴 NOT DONE | Q3–Q4 2026 | Diagnostics & validation |

### 🔴 **Phase 7–9: Complex Pathologies & ML** (PLANNED)
| Subtask | Status | Timeline |
|---------|--------|----------|
| **Phase 9: Global Circulation Model** | 🔴 NOT DONE | Q3 2026–Q2 2027 |
| — 3D/0D coupling | 🔴 NOT DONE | Q3 2026 |
| — Cholesterol plaque modeling | 🔴 NOT DONE | Q4 2026 |
| — Aortic rupture risk simulation | 🔴 NOT DONE | Q1 2027 |
| — Stent placement | 🔴 NOT DONE | Q2 2027 |
| — Hemolysis analysis | 🔴 NOT DONE | Q3 2027 |
| **Phase 10: Machine Learning** | 🟡 PARTIAL | Q1 2025–Q2 2027 |
| — Dataset (Kaggle) | ✅ DONE | Q1 2025 |
| — Health prediction model | 🔴 NOT DONE | Q1–Q4 2027 |
| — Diagnostic assistant (image recognition) | 🔴 NOT DONE | Q2–Q4 2027 |

---

## Critical Path & Dependencies

```
Phase 1–3 (Meshing) ✅
    ↓
Phase 4 (Aorta + CFD) ✅
    ↓
Phase 5 (0D Pathophysiology) ✅
    ├─→ Phase 7 (0D→2D Sync) ⚠️ BLOCKING ←─────┐
    │       (BC integration, WSS fields)           │
    │       Duration: 2–3 months                  │
    │       ↓                                      │
    └→ Phase 6 (FSI Prototype) 🔄 IN PROGRESS     │
            ├─→ Phase 8 (MFEM + FEBio) ⚠️ CRITICAL ├─ Path
            │       Benchmarks 1–3                │
            │       Turbine integration           │
            │       Advanced wall model           │
            │       Duration: 9–12 months         │
            │       ↓                              │
            └→ Phase 9 (Global Circulation) 🔴    │
                    (3D/0D coupling, pathology)   │
                    Duration: 6–9 months          │
                    ↓                              │
                    Phase 10 (ML models) 🔴       │
                            (Prediction, diagnosis)│
                            Duration: 6–12 months │
                            ↓
                    DELIVERY TARGET: Q4 2027 ←────┘
```

---

## Resource Budget & Constraints

### Smoke-Test Mode (Current)
- **Time step:** `DT = 2.5e-5 s`
- **Max steps:** `400`
- **Physical window:** `0.010 s` (≈ 1.2% of cardiac cycle)
- **CPU time target:** ≈ 15 minutes
- **Purpose:** Validation of BDF2 convergence and FSI residual stability

### Production Mode (Not Yet)
- **Full cardiac period:** `T = 0.8333 s`
- **Steps required:** 33,333 (at same DT)
- **Estimated CPU time:** ≈ 33 hours (single core)
- **Status:** ⏸️ Deferred until profiling is complete
- **Optimization plan:**
  1. Profile DolfinX / PETSc assembly, LU factorization, Newton iterations
  2. Improve algorithm without language change (preconditioning, MPI, adaptive I/O)
  3. Migrate to C++ (MFEM) with PETSc/Hypre for HPC
  4. Optional: GPU acceleration for hot kernels

---

## Technology Stack

| Layer | Tool | Purpose | Status |
|-------|------|---------|--------|
| **Meshing** | Gmsh | Geometry → mesh | ✅ Active |
| **CFD (2D/3D)** | DolfinX/FEniCSx | Navier–Stokes solver | ✅ Active |
| **Linear algebra** | PETSc | SNES, LU, preconditioning | ✅ Active |
| **0D Modeling** | Python (NumPy/SciPy) | Heart engine, Windkessel | ✅ Active |
| **FSI (next)** | FEBio | Hyperelastic tissue | 🔴 Not started |
| **FSI (production)** | MFEM | Monolithic coupled solver | 🔴 Not started |
| **ML (next)** | TensorFlow / PyTorch | Health prediction | 🔴 Not started |
| **Visualization** | ParaView + Matplotlib | Results / reports | ✅ Active |

---

## Key Milestones & Go/No-Go Decisions

### ✅ Completed Milestones
- **2024-06:** Mesher + geometry validation
- **2024-12:** IPCS solver stable, mass error < 1.5%
- **2025-06:** 5 WHO pathophysiology scenarios + 0D sweep complete

### 🟡 Upcoming Critical Decisions
1. **Q3 2025:** Phase 7 integration complete? → **GO** to FSI, **NO-GO** = resolve BC integration issues
2. **Q4 2025:** Benchmark 1 (2D channel FSI) validated? → **GO** to Benchmark 2, **NO-GO** = debug MFEM/FEBio setup
3. **Q2 2026:** Benchmark 3 (aortic arch FSI) stable? → **GO** to turbine, **NO-GO** = refine wall model
4. **Q3 2026:** Turbine + pulsatile flow working? → **GO** to pathology models, **NO-GO** = performance tuning
5. **Q4 2026:** Global circulation coupling stable? → **GO** to ML, **NO-GO** = resolve 3D/0D coupling
6. **Q2 2027:** ML models (prediction + diagnosis) trained and validated? → **DELIVERY READY**

---

## Risk Register

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| FSI convergence instability (Phase 8) | 🔴 High | Medium | Use smoke-test budget; pre-validate on 2D benchmarks |
| DolfinX performance scaling | 🔴 High | Medium | Profile early; migrate to MFEM if needed; HPC resources |
| FEBio integration complexity | 🔴 High | High | Start FEBio setup in Q3 2025; reference MLG_final_.pdf |
| Wall model (nonlinear, h(s)) too complex | 🟡 Medium | Medium | Use simplified model first; refine iteratively |
| ML dataset quality (Kaggle) | 🟡 Medium | Medium | Perform data audit early; consider synthetic data augmentation |
| Phase 7 BC integration delays | 🔴 High | High | **BLOCKING ITEM** — Assign dedicated resource NOW |

---

## Recommended Next Actions (Q4 2025)

1. **START NOW (October 2025):**
   - [ ] Assign dedicated developer to Phase 7 (0D→2D CFD sync)
   - [ ] Set up FEBio + MFEM development environment
   - [ ] Create detailed FSI solver specification document

2. **Q4 2025 Deliverables:**
   - [ ] Q_AO(t), P_AO(t) boundary conditions working in IPCS
   - [ ] WSS field comparison for ≥2 WHO scenarios (baseline + hypertension)
   - [ ] FEBio 2D channel mesh created and preprocessed

3. **Q1 2026 Target:**
   - [ ] Benchmark 1 (2D pulsating channel with elastic wall) running
   - [ ] MFEM/MFEMiFSI architecture document completed
   - [ ] Performance baseline: wall clock time for 100 time steps

---

**Last Updated:** October 2025  
**Project Owner:** Siemens Portfolio CFD & Healthcare  
**Questions?** Refer to `Plan.md` for detailed technical specs per phase.
