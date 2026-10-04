# Blood Pump CFD & Aortic FSI — Project Roadmap
**Internal R&D Project**  
**Project Start:** August 2026 | **Last Updated:** October 2026  
**Stack:** Python · DolfinX 0.11 · PETSc · Gmsh · ParaView

---

## Timeline (Gantt Chart)

```mermaid
gantt
    title Blood Pump CFD & Aortic FSI Simulation — Detailed Roadmap
    dateFormat YYYY-MM-DD
    axisFormat %b %Y
    
    section COMPLETED
    Phase 1: Geometry & Mesh          :done, p1, 2026-08-01, 31d
    Phase 2: Rigid CFD                :done, p2, 2026-08-15, 45d
    Phase 3: 0D Pathophysiology       :done, p3, 2026-09-01, 31d
    
    section IN PROGRESS
    Phase 4: FSI (Partitioned ALE)    :active, p4, 2026-09-15, 60d
    
    section PLANNED
    Phase 5: Pathological Runs        :crit, p5, 2026-11-01, 45d
    Phase 6: Stabilization            :crit, p6, 2026-12-15, 60d
    Phase 7: Visualization & Reports  :p7, 2027-02-01, 45d
    Phase 8: 3D Upgrade & Production  :crit, p8, 2027-03-15, 90d
    
    section Detailed Subtasks (Phase 4 — FSI)
    4.1: Submesh extraction           :done, sub4_1, 2026-09-15, 10d
    4.2: Interface DOF mapping (KDTree) :done, sub4_2, 2026-09-20, 15d
    4.3: Aitken relaxation loop       :done, sub4_3, 2026-09-30, 15d
    4.4: Zero-order predictor         :done, sub4_4, 2026-10-05, 10d
    4.5: VTK export & visualization   :done, sub4_5, 2026-10-10, 10d
    4.6: Detailed logging             :done, sub4_6, 2026-10-15, 5d
    
    section Detailed Subtasks (Phase 5 — Pathological Runs)
    5.1: Connect NPZ boundary conditions :crit, sub5_1, 2026-11-01, 15d
    5.2: Run 3 scenarios (baseline/hypertension/clot) :crit, sub5_2, 2026-11-15, 20d
    5.3: Extract WSS, pressure, displacement :sub5_3, 2026-12-01, 15d
    5.4: Build comparison plots       :sub5_4, 2026-12-10, 10d
    
    section Detailed Subtasks (Phase 6 — Stabilization)
    6.1: SUPG stabilization           :crit, sub6_1, 2026-12-15, 20d
    6.2: Nonlinear wall mechanics     :crit, sub6_2, 2027-01-05, 30d
    6.3: Smoother geometry (rounded)  :sub6_3, 2027-02-01, 20d
    6.4: Performance profiling        :sub6_4, 2027-02-15, 15d
    
    section Detailed Subtasks (Phase 7 — Visualization)
    7.1: ParaView animations          :sub7_1, 2027-02-01, 20d
    7.2: PDF report generation        :sub7_2, 2027-02-20, 15d
    7.3: Interactive HTML dashboard   :sub7_3, 2027-03-05, 15d
    7.4: Comparison tables & metrics  :sub7_4, 2027-03-15, 10d
    
    section Detailed Subtasks (Phase 8 — 3D & Production)
    8.1: 3D aorta geometry import      :crit, sub8_1, 2027-03-15, 30d
    8.2: MFEM / FEBio evaluation      :crit, sub8_2, 2027-04-15, 30d
    8.3: Monolithic FSI solver        :crit, sub8_3, 2027-05-15, 45d
    8.4: Large-scale runs (full cycle) :crit, sub8_4, 2027-07-01, 45d
    8.5: Validation & benchmarking    :sub8_5, 2027-08-15, 30d
    
    milestone M1: Phase 4 Complete (FSI stable), 2026-11-01, 1d
    milestone M2: Phase 5 Complete (3 scenarios), 2026-12-01, 1d
    milestone M3: Phase 6 Complete (SUPG + Neo-Hookean), 2027-02-15, 1d
    milestone M4: Phase 7 Complete (Visual reports), 2027-03-15, 1d
    milestone M5: 3D model & MFEM ready, 2027-06-01, 1d
    milestone M6: DELIVERY — Production FSI, 2027-09-01, 1d
```

---

## Phase Overview & Status

### ✅ **Phase 1: Geometry & Mesh** (COMPLETED — August 2026)

| Item | File | Status |
|------|------|--------|
| Parametric 2D aorta | `build_schematic_2d.py` | ✅ DONE |
| Base mesh (rigid) | `schematic_2d.msh` | ✅ DONE |
| Refined geometry | `schematic_rects.py` | ✅ DONE |
| FSI mesh | `schematic_rects_fsi.msh` | ✅ DONE |
| **Result** | 3683 nodes, 7792 elements, inlet + 3 outlets | ✅ |

---

### ✅ **Phase 2: Rigid CFD** (COMPLETED — August–September 2026)

| Task | Implementation | Status | Notes |
|------|-----------------|--------|-------|
| IPCS scheme | `ns_aorta_4outlet_heart_pulsatile.py` | ✅ DONE | Pulsatile flow, 4 outlets |
| RC Windkessel | `ns_aorta_4outlet_heart_pulsatile_RC_fixed.py` | ✅ DONE | Outlet impedance modeling |
| BDF2 scheme | `aorta_bdf2.py` | ✅ DONE | Final validated solver |
| 3-scheme comparison | `compare_aorta_3schemes.py` | ✅ DONE | IPCS vs Newton-CN vs BDF2 |
| 0D heart model | `heart_engine.py` | ✅ DONE | Time-varying elastance |
| **IPCS fixes applied** | `apply_lifting` + skew-symmetric convection | ✅ DONE | Mass error < 1.5% |

---

### ✅ **Phase 3: 0D Pathophysiology** (COMPLETED — September 2026)

| Scenario | File | Status | Output |
|----------|------|--------|--------|
| Baseline | `turbine_bc_0_baseline.npz` | ✅ DONE | Normal hemodynamics |
| Mild hypertension | `turbine_bc_1_hypertension_mild.npz` | ✅ DONE | Elevated R_sys |
| Hypertension | `turbine_bc_1_hypertension.npz` | ✅ DONE | High WSS trigger |
| Diabetes (stiff arteries) | `turbine_bc_2_diabetes_stiff.npz` | ✅ DONE | Reduced C_ao |
| Aortic stenosis | `turbine_bc_3_aortic_stenosis.npz` | ✅ DONE | High R_aortic |
| Cholesterol stenosis | `turbine_bc_3_cholesterol_stenosis.npz` | ✅ DONE | Plaque blockage |
| Severe combined | `turbine_bc_4_severe_combined.npz` | ✅ DONE | Multi-pathology |
| WHO killer combo | `turbine_bc_4_who_killer_combo.npz` | ✅ DONE | HTN + DM + plaque |
| Heart attack | `turbine_bc_5_heart_attack.npz` | ✅ DONE | Reduced ejection |
| Stroke | `turbine_bc_6_stroke.npz` | ✅ DONE | Clot + low flow |
| Cardiac arrest | `turbine_bc_7_cardiac_arrest.npz` | ✅ DONE | Minimal output |
| Clot (mild/mod/sev) | `turbine_bc_clot_*.npz` | ✅ DONE | Thrombus progression |
| Atrial fibrillation | `turbine_bc_clot_afib.npz` | ✅ DONE | Irregular rhythm |
| **Analysis** | `markov_cfd_analysis.py` | ✅ DONE | Progression modeling |
| **Dashboards** | `dashboard_*.pdf` | ✅ DONE | Visual reports |

---

### 🔄 **Phase 4: FSI (Partitioned ALE)** (IN PROGRESS — Sep–Oct 2026)

**Goal:** Strongly-coupled two-way FSI on elastic aortic wall.

#### Solver Architecture
```
┌─────────────────────────────────────────────┐
│  Fluid Domain (BDF2 ALE Navier-Stokes)      │
│  ∂u/∂t + (u·∇)u + ∇p = ν∇²u,  ∇·u = 0    │
└─────────────────────────────────────────────┘
         ↓ traction: σ·n
┌─────────────────────────────────────────────┐
│  Interface DOF Mapping (KDTree)             │
│  Error < 1e-9 m                             │
└─────────────────────────────────────────────┘
         ↓ displacement: d
┌─────────────────────────────────────────────┐
│  Solid Domain (Plane Strain Elasticity)     │
│  ρ_s ∂²d/∂t² = ∇·σ_s + f                   │
└─────────────────────────────────────────────┘
         ↓ 
┌─────────────────────────────────────────────┐
│  ALE Mesh Smoothing (Laplace)               │
│  ∇²X_ALE = 0  (boundary: d)                │
└─────────────────────────────────────────────┘
         ↓
┌─────────────────────────────────────────────┐
│  Aitken Relaxation (dynamic ω)              │
│  ω ∈ [0.03, 0.80]                           │
└─────────────────────────────────────────────┘
```

| Component | Implementation | Status |
|-----------|-----------------|--------|
| Fluid mesh extraction | `aorta_fsi_bdf2.py` | ✅ DONE |
| Solid mesh extraction | KDTree submesh | ✅ DONE |
| Interface mapping | K-d tree (error < 1e-9 m) | ✅ DONE |
| Aitken loop | Dynamic ω, convergence check | ✅ DONE |
| Zero-order predictor | Eliminates velocity shocks | ✅ DONE |
| VTK output | 40 frames per 400 steps | ✅ DONE |
| Detailed logging | Step, FSI-iter, traction, pressure, displacement | ✅ DONE |

#### Current Parameters
```
DT                 = 2.5e-4 s
MAX_STEPS          = 400
WRITE_EVERY        = 10
FSI_RELAX_MAX      = 0.80
MAX_INTERFACE_DISP = 2.0e-3 m (2 mm)
```

#### Known Limitations (Accepted for Phase 4)
- 2D Plane Strain: branches bend like cantilevers (not radial expansion)
- Numerical noise at high Re (needs SUPG in Phase 6)
- Sharp corners cause velocity singularities (smoothed in Phase 6)

---

### 🔴 **Phase 5: Pathological Runs & Analysis** (PLANNED — Nov–Dec 2026)

**Goal:** Run FSI for 3+ clinical scenarios, extract hemodynamic metrics.

| Task | Target | Timeline | Notes |
|------|--------|----------|-------|
| Connect NPZ boundary conditions | turbine_bc_*.npz → inlet in FSI solver | Nov 1–15 | Link Phase 3 to Phase 4 |
| Run baseline | Normal aorta FSI | Nov 15–30 | Reference case |
| Run hypertension | Elevated R_sys FSI | Nov 20–Dec 5 | High WSS |
| Run clot_severe | Thrombus blockage FSI | Nov 25–Dec 10 | Low WSS, stagnation |
| Export metrics | WSS, pressure, displacement to CSV | Dec 1–15 | Post-processing |
| Build comparison plots | Matplotlib figures (baseline vs pathologies) | Dec 10–20 | Visual analysis |
| Write clinical interpretation | Case-by-case summary | Dec 15–25 | Medical context |

---

### 🔴 **Phase 6: Stabilization & Quality** (PLANNED — Dec 2026 – Feb 2027)

**Goal:** Improve solver robustness, physics accuracy, geometry quality.

| Task | Technology | Timeline | Impact |
|------|-----------|----------|--------|
| SUPG stabilization | Streamline Upwind Petrov-Galerkin | Dec 15 – Jan 5 | Reduce numerical noise at high Re |
| Nonlinear wall mechanics | Neo-Hookean hyperelasticity | Jan 5 – Feb 5 | Replace linear elasticity |
| Smoother geometry | Rounded branch corners | Feb 1 – 15 | Eliminate velocity singularities |
| Performance profiling | PETSc assembly, LU timing | Feb 15 – Mar 1 | Baseline for HPC scaling |

---

### 🔴 **Phase 7: Visualization & Reporting** (PLANNED — Feb–Mar 2027)

**Goal:** Create publication-ready dashboards and reports.

| Deliverable | Tool | Timeline |
|-------------|------|----------|
| ParaView animations | `.avi` files per scenario | Feb 1–20 |
| PDF technical report | Scheme description + results | Feb 20 – Mar 5 |
| Interactive HTML dashboard | Web UI with plots/tables | Mar 5–15 |
| Comparison metrics table | WSS, pressure, compliance | Mar 10–20 |
| Presentation slides | Executive summary | Mar 15–25 |

---

### 🔴 **Phase 8: 3D Upgrade & Production** (PLANNED — Mar–Aug 2027)

**Goal:** Scale to 3D, migrate to production solver (MFEM/FEBio), full cardiac cycle.

| Task | Timeline | Deliverable |
|------|----------|-------------|
| 3D aorta geometry import (CT/MRI) | Mar 15 – Apr 15 | Realistic 3D model |
| MFEM / FEBio evaluation | Apr 15 – May 15 | Monolithic FSI framework |
| Monolithic FSI solver implementation | May 15 – Jun 30 | Production-ready code |
| Large-scale runs (full T=0.833 s cardiac cycle) | Jul 1 – Aug 15 | Complete hemodynamic cycle |
| Validation & benchmarking | Aug 15 – Sep 1 | Comparison vs in-vitro data |

---

## Critical Path & Dependencies

```
Phase 1 (Geometry) ✅
    ↓
Phase 2 (Rigid CFD) ✅
    ↓
Phase 3 (0D Pathophysiology) ✅
    ↓
Phase 4 (FSI — Partitioned) 🔄
    ├─→ Phase 5 (Pathological Runs) 🔴 ← BLOCKING
    │   (must connect NPZ → inlet)
    │   Duration: 6 weeks
    │   ↓
    └→ Phase 6 (Stabilization) 🔴
        (SUPG + Neo-Hookean + smoothing)
        Duration: 8 weeks
        ↓
    Phase 7 (Visualization) 🔴
    (PDF + HTML dashboards)
    Duration: 6 weeks
        ↓
    Phase 8 (3D + Production) 🔴
    (MFEM, full cycle)
    Duration: 22 weeks
        ↓
    🎯 DELIVERY: Q3 2027
```

---

## Key Milestones

| Milestone | Target Date | Criterion |
|-----------|-------------|-----------|
| **M1** Phase 4 stable FSI | Nov 1, 2026 | Zero-order predictor validated, Aitken ω ∈ [0.03, 0.80] |
| **M2** Phase 5 complete (3 scenarios) | Dec 1, 2026 | WSS + pressure + displacement exported, plots done |
| **M3** Phase 6 stabilization | Feb 15, 2027 | SUPG implemented, Neo-Hookean elastic wall works |
| **M4** Phase 7 reports | Mar 15, 2027 | Interactive HTML + PDF dashboards ready |
| **M5** 3D + MFEM ready | Jun 1, 2027 | 3D model imported, MFEM framework evaluated |
| **M6** 🎯 DELIVERY | Sep 1, 2027 | Production FSI, full cardiac cycle, validation complete |

---

## Technology Stack

| Layer | Tool | Version | Status |
|-------|------|---------|--------|
| **Mesh** | Gmsh | 4.13 | ✅ Active |
| **CFD/FSI** | DolfinX | 0.11 | ✅ Active |
| **Linear Algebra** | PETSc | 3.25 | ✅ Active |
| **0D Models** | Python/NumPy/SciPy | 3.11 | ✅ Active |
| **Visualization** | ParaView | 5.13 | ✅ Active |
| **Production FSI** | MFEM | next | 🔴 Phase 8 |
| **Tissue Mechanics** | FEBio | next | 🔴 Phase 8 |
| **ML/Analytics** | TensorFlow/PyTorch | future | 🔴 Beyond scope |

---

## Risk Register

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| FSI divergence on pathological cases (Phase 5) | 🔴 High | Medium | Pre-test with simple stenosis; use smaller time step |
| Aitken ω saturation at 0.80 (slow convergence) | 🟡 Medium | Medium | Consider Anderson mixing as alternative |
| SUPG implementation complexity (Phase 6) | 🟡 Medium | Low | Reference FEniCSx documentation; test on 2D cylinder |
| 3D mesh quality (Phase 8) | 🔴 High | Medium | Validate with industry CT scans early; use mesh quality metrics |
| MFEM monolithic solver instability | 🔴 High | Medium | Benchmark against DolfinX on 2D first; staged migration |

---

## Recommended Actions (October 2026)

### Immediate (This Week)
- [ ] Confirm Phase 4 FSI stability on baseline case
- [ ] Run 50-step test with clot_severe BC
- [ ] Document any convergence issues

### Next 2 Weeks
- [ ] Plan Phase 5 BC integration (connect turbine_bc_*.npz → inlet)
- [ ] Create test harness for multi-scenario runs
- [ ] Set up CSV export for WSS/pressure/displacement

### Next Month
- [ ] Complete Phase 5: Run all 3 scenarios (baseline, HTN, clot)
- [ ] Extract comparison plots
- [ ] Write preliminary clinical interpretation

### Strategic (Q1 2027)
- [ ] Begin Phase 6: SUPG stabilization + Neo-Hookean
- [ ] Plan Phase 8: 3D geometry acquisition
- [ ] Evaluate MFEM / FEBio licensing & integration

---

**Last Updated:** October 2026  
**Project Owner:** CFD / FSI Team  
**Status:** On track for Q3 2027 delivery
