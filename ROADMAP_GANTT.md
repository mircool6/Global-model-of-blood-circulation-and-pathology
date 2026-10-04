# Blood Pump CFD & Aortic FSI — Project Roadmap
**Internal R&D Project**  
**Project Start:** July 2026 | **Last Updated:** October 2026  
**Project Owner:** Artem Voitenko  
**Stack:** Python · DolfinX 0.11 · PETSc · Gmsh · ParaView · FreeCAD

---

## Timeline (Gantt Chart)

```mermaid
gantt
    title Blood Pump CFD & Aortic FSI Simulation — Detailed Roadmap
    dateFormat YYYY-MM-DD
    axisFormat %b %Y

    section COMPLETED
    Phase 0: Engine Development       :done, p0, 2026-07-01, 45d
    Phase 1: Geometry & Mesh          :done, p1, 2026-08-15, 31d
    Phase 2: Rigid CFD                :done, p2, 2026-09-01, 45d
    Phase 3: Data Parsing & Analytics :done, p3, 2026-09-15, 45d

    section IN PROGRESS
    Phase 4: FSI (Partitioned ALE)    :active, p4, 2026-10-01, 60d

    section PLANNED
    Phase 5: Pathological Runs        :crit, p5, 2026-11-15, 45d
    Phase 6: Stabilization            :crit, p6, 2026-12-30, 60d
    Phase 7: Visualization & Reports  :p7, 2027-02-15, 45d
    Phase 8: 3D Upgrade & Production  :crit, p8, 2027-03-30, 90d

    section Detailed Subtasks (Phase 0 — Engine Development)
    0.1: FreeCAD parameterization     :done, sub0_1, 2026-07-01, 15d
    0.2: Case generation engine       :done, sub0_2, 2026-07-15, 20d
    0.3: 9 simulation scenarios        :done, sub0_3, 2026-08-01, 10d

    section Detailed Subtasks (Phase 1 — Geometry)
    1.1: Base aorta scheme            :done, sub1_1, 2026-08-15, 15d
    1.2: Refined mesh generation      :done, sub1_2, 2026-08-25, 10d
    1.3: FSI mesh creation            :done, sub1_3, 2026-09-01, 6d

    section Detailed Subtasks (Phase 2 — Rigid CFD)
    2.1: IPCS / BDF2 solver           :done, sub2_1, 2026-09-01, 20d
    2.2: Windkessel BCs               :done, sub2_2, 2026-09-15, 15d
    2.3: Mass-error corrections       :done, sub2_3, 2026-09-25, 10d

    section Detailed Subtasks (Phase 3 — Data Parsing & Analytics)
    3.1: WHO / GBD / PubMed parser   :done, sub3_1, 2026-09-15, 15d
    3.2: 20 NPZ scenario exports      :done, sub3_2, 2026-09-25, 10d
    3.3: Dashboard & PDF reports      :done, sub3_3, 2026-10-05, 10d

    section Detailed Subtasks (Phase 4 — FSI)
    4.1: Submesh extraction           :done, sub4_1, 2026-10-01, 10d
    4.2: Interface DOF mapping        :done, sub4_2, 2026-10-08, 15d
    4.3: Aitken relaxation loop       :done, sub4_3, 2026-10-18, 15d
    4.4: Zero-order predictor         :done, sub4_4, 2026-10-28, 10d
    4.5: VTK export & visualization   :done, sub4_5, 2026-11-01, 10d
    4.6: Detailed logging             :done, sub4_6, 2026-11-08, 5d

    section Detailed Subtasks (Phase 5 — Pathological Runs)
    5.1: Connect NPZ boundary conditions :crit, sub5_1, 2026-11-15, 15d
    5.2: Run 3 scenarios (baseline/hypertension/clot) :crit, sub5_2, 2026-11-30, 20d
    5.3: Extract WSS, pressure, displacement :sub5_3, 2026-12-15, 15d
    5.4: Build comparison plots       :sub5_4, 2026-12-25, 10d

    section Detailed Subtasks (Phase 6 — Stabilization)
    6.1: SUPG stabilization           :crit, sub6_1, 2026-12-30, 20d
    6.2: Nonlinear wall mechanics     :crit, sub6_2, 2027-01-20, 30d
    6.3: Smoother geometry (rounded)  :sub6_3, 2027-02-15, 20d
    6.4: Performance profiling        :sub6_4, 2027-03-01, 15d

    section Detailed Subtasks (Phase 7 — Visualization)
    7.1: ParaView animations          :sub7_1, 2027-02-15, 20d
    7.2: PDF report generation        :sub7_2, 2027-03-05, 15d
    7.3: Interactive HTML dashboard   :sub7_3, 2027-03-20, 15d
    7.4: Comparison tables & metrics  :sub7_4, 2027-04-01, 10d

    section Detailed Subtasks (Phase 8 — 3D & Production)
    8.1: 3D geometry import           :crit, sub8_1, 2027-03-30, 30d
    8.2: MFEM / FEBio evaluation      :crit, sub8_2, 2027-04-30, 30d
    8.3: Monolithic FSI solver        :crit, sub8_3, 2027-05-30, 45d
    8.4: Large-scale runs (full cycle) :crit, sub8_4, 2027-07-15, 45d
    8.5: Validation & benchmarking    :sub8_5, 2027-08-30, 30d

    milestone M0: Engine setup validated, 2026-08-15, 1d
    milestone M1: Phase 4 Complete (FSI stable), 2026-12-01, 1d
    milestone M2: Phase 5 Complete (3 scenarios), 2027-01-01, 1d
    milestone M3: Phase 6 Complete (SUPG + Neo-Hookean), 2027-03-15, 1d
    milestone M4: Phase 7 Complete (Visual reports), 2027-04-15, 1d
    milestone M5: 3D model & MFEM ready, 2027-06-15, 1d
    milestone M6: DELIVERY — Production FSI, 2027-09-01, 1d
```

---

## Phase Overview & Status

### ✅ **Phase 0: Engine Development** (COMPLETED — July 2026)

| Item | Scope | Status |
|------|-------|--------|
| FreeCAD-based geometry engine | Parameterized geometry generation | ✅ DONE |
| Simulation case generator | Multi-case scenario workflow | ✅ DONE |
| Scenario family | 9 core cases for engine / pressure / flow studies | ✅ DONE |
| **Result** | Automated setup for systematic CFD and FSI case generation | ✅ |

---

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

### ✅ **Phase 3: Data Parsing & Clinical Analytics** (COMPLETED — September 2026)

| Component | Scope | Status | Output |
|-----------|-------|--------|--------|
| WHO / GBD / PubMed parser | Clinical and risk-factor data ingestion | ✅ DONE | Structured dataset |
| Scenario export | 20 NPZ boundary-condition files | ✅ DONE | `turbine_bc_*.npz` |
| Dashboard generation | PDF, reports, comparative plots | ✅ DONE | Visualization outputs |
| **Result** | Clinical parameterization for simulation scenarios | ✅ |

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
| Connect NPZ boundary conditions | turbine_bc_*.npz → inlet in FSI solver | Nov 15–30 | Link Phase 3 to Phase 4 |
| Run baseline | Normal aorta FSI | Nov 30–Dec 10 | Reference case |
| Run hypertension | Elevated R_sys FSI | Dec 1–15 | High WSS |
| Run clot_severe | Thrombus blockage FSI | Dec 10–20 | Low WSS, stagnation |
| Export metrics | WSS, pressure, displacement to CSV | Dec 15–30 | Post-processing |
| Build comparison plots | Matplotlib figures (baseline vs pathologies) | Dec 20–Jan 5 | Visual analysis |
| Write clinical interpretation | Case-by-case summary | Jan 5–15 | Medical context |

---

### 🔴 **Phase 6: Stabilization & Quality** (PLANNED — Dec 2026 – Feb 2027)

**Goal:** Improve solver robustness, physics accuracy, geometry quality.

| Task | Technology | Timeline | Impact |
|------|-----------|----------|--------|
| SUPG stabilization | Streamline Upwind Petrov-Galerkin | Dec 30 – Jan 20 | Reduce numerical noise at high Re |
| Nonlinear wall mechanics | Neo-Hookean hyperelasticity | Jan 20 – Feb 20 | Replace linear elasticity |
| Smoother geometry | Rounded branch corners | Feb 15 – Mar 1 | Eliminate velocity singularities |
| Performance profiling | PETSc assembly, LU timing | Mar 1 – 15 | Baseline for HPC scaling |

---

### 🔴 **Phase 7: Visualization & Reporting** (PLANNED — Feb–Mar 2027)

**Goal:** Create publication-ready dashboards and reports.

| Deliverable | Tool | Timeline |
|-------------|------|----------|
| ParaView animations | `.avi` files per scenario | Feb 15–Mar 5 |
| PDF technical report | Scheme description + results | Mar 5–20 |
| Interactive HTML dashboard | Web UI with plots/tables | Mar 20–Apr 1 |
| Comparison metrics table | WSS, pressure, compliance | Apr 1–10 |
| Presentation slides | Executive summary | Apr 10–20 |

---

### 🔴 **Phase 8: 3D Upgrade & Production** (PLANNED — Mar–Aug 2027)

**Goal:** Scale to 3D, migrate to production solver (MFEM/FEBio), full cardiac cycle.

| Task | Timeline | Deliverable |
|------|----------|-------------|
| 3D aorta geometry import (CT/MRI) | Mar 30 – Apr 30 | Realistic 3D model |
| MFEM / FEBio evaluation | Apr 30 – May 30 | Monolithic FSI framework |
| Monolithic FSI solver implementation | May 30 – Jul 15 | Production-ready code |
| Large-scale runs (full T=0.833 s cardiac cycle) | Jul 15 – Aug 30 | Complete hemodynamic cycle |
| Validation & benchmarking | Aug 30 – Sep 15 | Comparison vs in-vitro data |

---

## Critical Path & Dependencies

```
Phase 0 (Engine Development) ✅
    ↓
Phase 1 (Geometry) ✅
    ↓
Phase 2 (Rigid CFD) ✅
    ↓
Phase 3 (Data Parsing & Analytics) ✅
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
| **M0** Engine setup validated | Aug 15, 2026 | FreeCAD engine + scenario workflow ready |
| **M1** Phase 4 stable FSI | Dec 1, 2026 | Zero-order predictor validated, Aitken ω ∈ [0.03, 0.80] |
| **M2** Phase 5 complete (3 scenarios) | Jan 1, 2027 | WSS + pressure + displacement exported, plots done |
| **M3** Phase 6 stabilization | Mar 15, 2027 | SUPG implemented, Neo-Hookean elastic wall works |
| **M4** Phase 7 reports | Apr 15, 2027 | Interactive HTML + PDF dashboards ready |
| **M5** 3D + MFEM ready | Jun 15, 2027 | 3D model imported, MFEM framework evaluated |
| **M6** 🎯 DELIVERY | Sep 1, 2027 | Production FSI, full cardiac cycle, validation complete |

---

## Technology Stack

| Layer | Tool | Version | Status |
|-------|------|---------|--------|
| **Geometry / Parameterization** | FreeCAD | 1.0 | ✅ Active |
| **Mesh** | Gmsh | 4.13 | ✅ Active |
| **CFD/FSI** | DolfinX | 0.11 | ✅ Active |
| **Linear Algebra** | PETSc | 3.25 | ✅ Active |
| **0D Models** | Python/NumPy/SciPy | 3.11 | ✅ Active |
| **Data Processing** | Pandas / NumPy / parsing scripts | current | ✅ Active |
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

**Status:** On track for Q3 2027 delivery

