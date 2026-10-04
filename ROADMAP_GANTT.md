# Blood Pump CFD & Aortic FSI — Project Roadmap
**Internal R&D Project**  
**Project Start:** July 2026 | **Last Updated:** October 2026  
**Project Owner:** Artem Voitenko  
**Stack:** Python · DolfinX 0.11 · PETSc · Gmsh · ParaView · FreeCAD · NumPy · SciPy · Pandas

---

## Executive Summary

This project develops a **two-way fluid-structure interaction (FSI) simulator** for hemodynamic modeling of aortic blood flow under normal and pathological conditions. The workflow progresses from 2D parametric geometry → rigid CFD validation → partitioned FSI coupling → clinical scenario analysis → 3D production system with full cardiac cycle simulation.

**Key Deliverables:**
- ✅ Partitioned FSI solver with Aitken acceleration (Phase 4)
- 🔴 Clinical scenario suite: baseline, hypertension, thrombosis (Phase 5)
- 🔴 Production-grade 3D solver with MFEM/FEBio (Phase 8)

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

| Item | Scope | Status | Key Files |
|------|-------|--------|-----------|
| FreeCAD-based geometry engine | Parameterized geometry generation | ✅ DONE | `freecad_*.py` |
| Simulation case generator | Multi-case scenario workflow | ✅ DONE | `case_generator.py` |
| Scenario family | 9 core cases for engine / pressure / flow studies | ✅ DONE | 9 `.FCStd` configs |
| **Result** | Automated setup for systematic CFD and FSI case generation | ✅ | — |

**Deliverables:**
- ✅ Parameterized FreeCAD models (9 geometry variants)
- ✅ Case setup automation (inlet pressure, outlet resistance swept)
- ✅ Scenario export (baseline, ±10% variations)

---

### ✅ **Phase 1: Geometry & Mesh** (COMPLETED — August 2026)

| Item | File | Spec | Status |
|------|------|------|--------|
| Parametric 2D aorta | `build_schematic_2d.py` | Inlet + 3 branches | ✅ DONE |
| Base mesh (rigid) | `schematic_2d.msh` | Gmsh 4.13 format | ✅ DONE |
| Refined geometry | `schematic_rects.py` | Edge refinement | ✅ DONE |
| FSI mesh | `schematic_rects_fsi.msh` | Fluid + solid partitions | ✅ DONE |
| **Result** | 3683 nodes, 7792 elements | Inlet + 3 outlets | ✅ |

**Mesh Details:**
```
Domain: 2D plane strain, aorta-like branching
Fluid mesh:    ~3200 elements, Re ≈ 600 at peak systole
Solid mesh:    ~500 elements, wall thickness ≈ 1 mm (scaled)
Interface:     ~80 DOF, smooth radial boundary
```

---

### ✅ **Phase 2: Rigid CFD** (COMPLETED — August–September 2026)

| Task | Implementation | Status | Scheme | Notes |
|------|-----------------|--------|--------|-------|
| IPCS scheme | `ns_aorta_4outlet_heart_pulsatile.py` | ✅ DONE | Fractional step | Pulsatile, 4 outlets |
| RC Windkessel | `ns_aorta_4outlet_heart_pulsatile_RC_fixed.py` | ✅ DONE | Outlet BC | Impedance R_prox, R_dist, C |
| BDF2 scheme | `aorta_bdf2.py` | ✅ DONE | 2nd-order time | Final validated solver |
| 3-scheme comparison | `compare_aorta_3schemes.py` | ✅ DONE | Benchmark | IPCS vs Newton-CN vs BDF2 |
| 0D heart model | `heart_engine.py` | ✅ DONE | Elastance | Time-varying elastance (0D) |
| **IPCS fixes applied** | `apply_lifting` + skew-symmetric convection | ✅ DONE | Numerical | Mass error < 1.5% |

**Solver Parameters:**
```
Time:      T = 0.833 s (1 cardiac cycle, 80 bpm)
Time step: Δt = 2e-4 s, 4165 steps per cycle
Pressure:  IPCS = Independent PP + Consistent Scaling
Velocity:  BDF2 = 2nd-order backward differentiation
Outlet BC: 4-element RC network per branch (6 unknowns)
```

**Validation:**
- ✅ Mass balance < 1.5% error
- ✅ Peak systolic flow ≈ 25 mL/s (physiologic)
- ✅ Diastolic pressure decay (Windkessel response)

---

### ✅ **Phase 3: Data Parsing & Clinical Analytics** (COMPLETED — September 2026)

| Component | Scope | Status | Output |
|-----------|-------|--------|--------|
| WHO / GBD / PubMed parser | Clinical & risk-factor data | ✅ DONE | Structured JSON |
| Scenario export | 20 NPZ boundary-condition files | ✅ DONE | `turbine_bc_*.npz` |
| Dashboard generation | PDF, reports, comparative plots | ✅ DONE | Matplotlib figures |
| **Result** | Clinical parameterization for simulation | ✅ | — |

**Dataset:**
```
Baseline:       Normal aortic properties, systolic 120 mmHg
Hypertension:   Elevated systolic (160 mmHg), increased R_sys
Thrombosis:     75% lumen blockage, stenosis at branch junction
Atherosclerosis: Diffuse wall stiffening (E_wall +30%)
Aneurysm:       Localized wall dilation (diameter +40%)

Export: 20 × 1D boundary conditions (inlet Q(t), outlet R, C, L)
        stored as NPZ for fast loading in Phase 5
```

**Key Metrics Tracked:**
- Inlet pulsatile flow rate Q(t)
- Outlet impedances (resistance, compliance, inertance)
- Wall elasticity parameters
- Risk factors (age, BMI, smoking status)

---

### 🔄 **Phase 4: FSI (Partitioned ALE)** (IN PROGRESS — Oct–Nov 2026)

**Goal:** Strongly-coupled two-way FSI on elastic aortic wall using partitioned approach with dynamic Aitken acceleration.

#### Solver Architecture
```
┌─────────────────────────────────────────────┐
│  Fluid Domain (BDF2 ALE Navier-Stokes)      │
│  ∂u/∂t + (u·∇)u + ∇p = ν∇²u,  ∇·u = 0    │
│  Pulsatile inlet + 4-outlet Windkessel BC   │
└─────────────────────────────────────────────┘
         ↓ traction: σ·n (stress → wall)
┌─────────────────────────────────────────────┐
│  Interface DOF Mapping (KDTree)             │
│  Error < 1e-9 m, projection by proximity    │
└─────────────────────────────────────────────┘
         ↓ displacement: d (wall → mesh)
┌─────────────────────────────────────────────┐
│  Solid Domain (Plane Strain Elasticity)     │
│  ρ_s ∂²d/∂t² = ∇·σ_s + f                   │
│  Linear: E_wall = 0.5 MPa, ν_wall = 0.3    │
└─────────────────────────────────────────────┘
         ↓ ALE mesh update
┌─────────────────────────────────────────────┐
│  ALE Mesh Smoothing (Laplace)               │
│  ∇²X_ALE = 0  (boundary: d)                │
│  Ensures element quality in fluid domain    │
└─────────────────────────────────────────────┘
         ↓
┌─────────────────────────────────────────────┐
│  Aitken Relaxation (dynamic ω)              │
│  ω ∈ [0.03, 0.80], convergence check       │
│  Typically 5–8 iterations per time step     │
└─────────────────────────────────────────────┘
```

#### Implementation Details

| Component | Implementation | Status | Details |
|-----------|-----------------|--------|---------|
| Fluid mesh extraction | `aorta_fsi_bdf2.py` | ✅ DONE | DolfinX MixedElement |
| Solid mesh extraction | KDTree submesh | ✅ DONE | Submesh from parent `msh` |
| Interface mapping | K-d tree (error < 1e-9 m) | ✅ DONE | Nearest-neighbor projection |
| Aitken loop | Dynamic ω, convergence check | ✅ DONE | Relative traction residual |
| Zero-order predictor | Eliminates velocity shocks | ✅ DONE | u_n extrapolation on mesh |
| VTK output | 40 frames per 400 steps | ✅ DONE | ParaView-compatible `.vtu` |
| Detailed logging | Step, FSI-iter, traction, pressure, displacement | ✅ DONE | CSV + terminal monitor |

#### Current Parameters
```
Time Integration:
  DT                 = 2.5e-4 s
  MAX_STEPS          = 400  (simulates ~0.1 s)
  WRITE_EVERY        = 10   (output every 2.5 ms)

FSI Coupling:
  FSI_RELAX_MAX      = 0.80   (Aitken acceleration factor)
  MAX_INTERFACE_DISP = 2.0e-3 m (2 mm, convergence threshold)
  CONVERGE_TOL       = 1e-5 (relative residual)

Physics:
  ν_fluid            = 3.5e-6 m²/s (blood kinematic viscosity)
  ρ_fluid            = 1050 kg/m³
  ρ_wall             = 1100 kg/m³
  E_wall             = 0.5 MPa (linear elasticity)
  ν_wall             = 0.3 (Poisson ratio)
```

#### Known Limitations (Accepted for Phase 4)
- ⚠️ **2D Plane Strain:** Branches bend like cantilevers (not true radial expansion)
  - *Workaround:* 3D upgrade in Phase 8
- ⚠️ **Numerical noise at high Re:** Needs upwind stabilization (SUPG in Phase 6)
  - *Impact:* ~±5% oscillation in wall shear stress near peak systole
- ⚠️ **Sharp corners:** Cause velocity singularities at branch junctions
  - *Mitigation:* Smooth geometry in Phase 6.3

#### Phase 4 Checkpoints
- [ ] ✅ **4.1–4.6:** Core FSI pipeline validated on baseline case
- [ ] **Nov 1:** VTK visualization of full 400-step run
- [ ] **Nov 8:** Convergence analysis (Aitken ω vs iteration count)
- [ ] **Next:** Ready for Phase 5 pathological cases

---

### 🔴 **Phase 5: Pathological Runs & Analysis** (PLANNED — Nov–Dec 2026)

**Goal:** Run FSI for 3+ clinical scenarios, extract hemodynamic metrics, validate solver on diverse physics.

**Scenarios:**
1. **Baseline (control):** Normal aortic properties, normal flow
2. **Hypertension:** Elevated outlet resistance, thick wall (E +20%)
3. **Thrombosis/Stenosis:** 75% lumen blockage, local pressure drop

#### Detailed Task Breakdown

| Task | Target | Timeline | Effort | Notes |
|------|--------|----------|--------|-------|
| **5.1** Connect NPZ boundary conditions | `turbine_bc_*.npz` → inlet in FSI solver | Nov 15–30 | 2 weeks | Link Phase 3 data to Phase 4 solver. Create adapter layer. Test on 1 case. |
| **5.2a** Run baseline FSI | 3 cardiac cycles (~2.5 s simulation) | Nov 30–Dec 10 | 1.5 weeks | Reference case. Debug logging. Check mass balance. |
| **5.2b** Run hypertension FSI | Elevated R_sys (×1.5), wall E (+20%) | Dec 1–15 | 2 weeks | High WSS regime. May need smaller Δt if unstable. |
| **5.2c** Run stenosis FSI | 75% blockage at branch 2, turbulent BC | Dec 10–20 | 2 weeks | Low WSS behind clot, high upstream. Risk of divergence. |
| **5.3** Extract hemodynamic metrics | WSS (Pa), pressure (mmHg), wall displacement (mm) | Dec 15–30 | 1.5 weeks | Post-process VTK. Export CSV per scenario. |
| **5.4a** Build comparison plots | Matplotlib: baseline vs HTN vs stenosis | Dec 20–Jan 5 | 1.5 weeks | WSS maps, pressure time-series, flow-rate overlay |
| **5.4b** Write clinical interpretation | Per-scenario summary, remarks on physics | Jan 5–15 | 1 week | Is WSS elevation in HTN clinically realistic? Discuss limitations. |

**Outputs:**
```
Phase 5 Deliverables:
├── NPZ adapters (3.py files)
├── 3 × FSI simulation runs (VTK archives)
├── CSV exports:
│   ├── baseline_WSS.csv, baseline_pressure.csv
│   ├── hypertension_*.csv
│   └── stenosis_*.csv
├── Matplotlib figures:
│   ├── WSS_map_comparison.png
│   ├── pressure_timeseries.png
│   └── flow_rate_overlay.png
└── Clinical_Report.md (1500 words)
```

**Risk:** FSI solver may diverge on stenosis case (high Re, recirculation). Mitigation: pre-test with Δt = 1e-4, use smaller Aitken ω.

---

### 🔴 **Phase 6: Stabilization & Quality** (PLANNED — Dec 2026 – Mar 2027)

**Goal:** Improve solver robustness, physics accuracy, and geometry quality to enable confident production runs.

| Task | Technology | Timeline | Impact | Complexity |
|------|-----------|----------|--------|------------|
| **6.1** SUPG stabilization | Streamline Upwind Petrov-Galerkin | Dec 30 – Jan 20 | Reduce numerical noise at high Re (~±2% vs ±5%) | Medium |
| **6.2** Nonlinear wall mechanics | Neo-Hookean hyperelasticity | Jan 20 – Feb 20 | Replace linear elasticity, handle large strains | High |
| **6.3** Smoother geometry | Rounded branch corners (fillets) | Feb 15 – Mar 1 | Eliminate velocity singularities | Low |
| **6.4** Performance profiling | PETSc assembly, LU timing, memory | Mar 1 – 15 | Baseline for HPC scaling (Phase 8) | Medium |

**6.1 SUPG Details:**
```
Goal: Add streamwise upwind diffusion to momentum equation
Implementation:
  1. Compute streamline direction: s = u / |u|
  2. Add artificial diffusion: τ_SUPG = 0.5 * h / (2|u|)
  3. Modify variational form: ∫ τ_SUPG (u·∇u) · v dx
Library support: FEniCSx documentation (form language)
Expected benefit: Reduce overshoots near high-gradient zones
```

**6.2 Neo-Hookean Details:**
```
Goal: Replace linear elasticity with hyperelasticity for realistic large-strain response
Strain energy: W = (μ/2)(tr(F^T F) - 3) - μ ln(det F) + (λ/2) ln²(det F)
  where F = I + ∇d (deformation gradient)
Implementation: ufl.ln(), ufl.sqrt() for symbolic derivatives
Wall parameters (literature values):
  μ ≈ 0.15 MPa (shear modulus)
  λ ≈ 0.5 MPa (Lamé parameter)
Expected behavior: Wall compliance ↓ at high pressure, mimics arterial stiffening
```

**6.3 Geometry Improvements:**
```
Current issues:
  ├─ Sharp 90° corners at branch junctions
  └─ Velocity singularities in CFD (u → ∞ at corners)

Fixes:
  ├─ Fillet radius r = 0.5 mm at all junctions
  ├─ Smooth inlet/outlet transitions
  └─ Re-mesh with Gmsh (refine near fillets)

Impact: Better convergence, realistic secondary flows
```

**6.4 Performance Baseline:**
```
Target metrics (per time step):
  ├─ Matrix assembly:        ~ 0.1 s
  ├─ LU factorization:       ~ 0.5 s
  ├─ Fluid solve:            ~ 1.5 s
  ├─ Solid solve:            ~ 0.3 s
  ├─ ALE update:             ~ 0.2 s
  └─ Total per step:         ~ 3.0 s (400 steps → 20 min wall clock)

Scaling targets for 3D (Phase 8):
  ├─ 2D:    ~3200 fluid DOF, ~500 solid DOF
  └─ 3D: ~100K fluid DOF, ~15K solid DOF (target: <50 s/step on 4-core)
```

---

### 🔴 **Phase 7: Visualization & Reporting** (PLANNED — Feb–Apr 2027)

**Goal:** Create publication-ready dashboards, videos, and technical reports.

| Deliverable | Tool | Format | Timeline | Audience |
|-------------|------|--------|----------|----------|
| **7.1** ParaView animations | `.mp4` (velocity, pressure, WSS) | Video per scenario | Feb 15–Mar 5 | Clinicians, reviewers |
| **7.2** PDF technical report | LaTeX + figures | Peer-review format | Mar 5–20 | Journal submission |
| **7.3** Interactive HTML dashboard | Plotly/Dash + JavaScript | Web UI | Mar 20–Apr 1 | Internal + web demo |
| **7.4** Comparison metrics table | Markdown + CSV | Static table | Apr 1–10 | Quick reference |
| **7.5** Presentation slides | PowerPoint / Beamer | Executive summary | Apr 10–20 | Stakeholders |

**7.1 ParaView Videos** (3 × scenario):
```
Each video (~3 min at 24 fps):
├─ Velocity contours + streamlines
├─ Pressure field + legend
├─ Wall shear stress (arrows on wall)
├─ Wall displacement magnitude
└─ Time stamp + cycle indicator
```

**7.2 PDF Report** (15–20 pages):
```
1. Introduction (clinical motivation)
2. Methods
   ├─ Geometry and meshing
   ├─ Equations (Navier-Stokes + elasticity)
   ├─ FSI coupling (partitioned ALE + Aitken)
   └─ Numerical parameters
3. Validation (rigid CFD vs literature)
4. Results
   ├─ Baseline case (WSS, pressure, compliance)
   ├─ Hypertension (comparison)
   └─ Stenosis (flow separation analysis)
5. Discussion
6. Limitations & future work
```

**7.3 HTML Dashboard:**
```
Components:
├─ Scenario selector (radio buttons)
├─ Time slider (0–0.833 s)
├─ Plot 1: Pressure vs time (Plotly.js)
├─ Plot 2: WSS heatmap (Plotly surface)
├─ Plot 3: Flow-rate vs scenario (Plotly bar)
├─ Download button (CSV, VTK)
└─ Embedded ParaView WebGL viewer (optional)
```

---

### 🔴 **Phase 8: 3D Upgrade & Production** (PLANNED — Mar–Aug 2027)

**Goal:** Scale solver to 3D, migrate to production framework (MFEM or FEBio), and run full cardiac cycle with validation.

| Task | Timeline | Technology | Deliverable | Complexity |
|------|----------|-----------|-------------|------------|
| **8.1** 3D geometry import | Mar 30 – Apr 30 | CT/MRI segmentation + CAD repair | 3D `.msh` file | High |
| **8.2** MFEM / FEBio evaluation | Apr 30 – May 30 | Proof-of-concept | Benchmark report | Medium |
| **8.3** Monolithic FSI solver | May 30 – Jul 15 | Coupled Newton solver | Production code | Very High |
| **8.4** Large-scale runs | Jul 15 – Aug 30 | Full T=0.833 s cycle (4000 steps) | 3D simulation archive | Medium |
| **8.5** Validation | Aug 30 – Sep 15 | Compare vs in-vitro PIV/MRI data | Technical report | Medium |

**8.1 3D Geometry Details:**
```
Data source:
  Option A: Open-access patient CT (public datasets)
  Option B: Synthetic 3D aorta (parametric CAD)

Workflow:
  1. Segment DICOM images (3D-Slicer or SimpleITK)
  2. Clean surface (remove artifacts, smooth)
  3. Export as STL
  4. Mesh with Gmsh (adaptive refinement)
  5. Create fluid/solid submeshes

Expected specs:
  ├─ Domain: ~4 cm aorta (ascending → thoracic)
  ├─ Fluid mesh: 150K–300K elements
  ├─ Solid mesh: 30K–50K elements (wall layer)
  └─ Simulation: ~4000 time steps (full cycle)
```

**8.2 MFEM vs FEBio Comparison:**
```
Criteria                 DolfinX (current)    MFEM                FEBio
──────────────────────────────────────────────────────────────────────
Monolithic FSI?          ❌ Partitioned        ✅ Yes             ✅ Yes
Open-source?             ✅ Yes               ✅ Yes             ⚠️ Free but proprietary
Nonlinear solvers?       ✅ SNESSolver        ✅ Yes             ✅ Yes
GPU support?             ⚠️ Limited          ✅ Yes             ⚠️ Limited
Learning curve?          Medium              High               Medium
```

**Recommendation:** MFEM for HPC scaling (GPU), FEBio as backup (commercial support).

**8.3 Monolithic FSI Solver:**
```
System:
  ┌─────────────────────┐
  │ Coupled residual:   │
  │ R_u = ... (NS)      │
  │ R_d = ... (elasticity)
  │ R_p = ... (pressure)
  └─────────────────────┘

Solve via:
  Newton-Raphson: x_{n+1} = x_n - J^{-1} R(x_n)
  where J = ∂R/∂x is Jacobian (3×3 block)
  
Advantages over partitioned:
  ✓ Better stability on stiff physics
  ✓ Fewer iterations per time step (~2–3 vs ~6–8)
  ✓ Scales better to 3D

Challenges:
  ✗ Larger linear system (N_fluid + N_solid)
  ✗ Jacobian more expensive to assemble
  ✗ Need good preconditioner
```

**8.4 Full-Cycle Simulation:**
```
Goal: Run complete 0.833 s cardiac cycle (80 bpm)

Timeline:
  ├─ Systole (0–0.3 s):   Rapid pressure rise, high flow
  ├─ Diastole (0.3–0.8 s): Gradual decay, backflow from outlets
  └─ Final (0.8–0.833 s):  Atrial refill

Outputs:
  ├─ VTK frames (every 10 ms)
  ├─ Time-series CSV (WSS, Q, P per node/segment)
  ├─ Compliance curves (ΔV/ΔP)
  └─ Phase-plane diagrams (flow vs pressure)
```

**8.5 Validation Against Experiments:**
```
Data sources:
  ├─ In-vitro PIV (Particle Image Velocimetry) — velocity maps
  ├─ In-vivo MRI (4D flow MRI) — patient data
  └─ Arterial compliance measurements — pressure-volume curves

Metrics:
  ├─ MAE(velocity):     < 10% of peak
  ├─ MAE(pressure):     < 5 mmHg
  ├─ Compliance error:  < 15%
  └─ WSS correlation:   R > 0.85
```

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
Phase 4 (FSI — Partitioned) 🔄 ← IN PROGRESS
    │
    ├─→ Phase 5 (Pathological Runs) 🔴
    │   Duration: 6 weeks
    │   Blocking: Phase 5.1 (BC integration)
    │   ↓
    │   Phase 6 (Stabilization) 🔴
    │   Duration: 8 weeks
    │   ↓
    │   Phase 7 (Visualization) 🔴
    │   Duration: 6 weeks
    │   ↓
    │   Phase 8 (3D + Production) 🔴
    │   Duration: 22 weeks
    │
    └──→ 🎯 DELIVERY: Sep 1, 2027
```

**Critical dependencies:**
- Phase 5 **blocks** Phase 6 (need validated pathological solutions first)
- Phase 6.1 (SUPG) **enables** confident Phase 8 3D runs
- Phase 8.3 (monolithic FSI) **requires** Phase 6.2 (Neo-Hookean elasticity)

---

## Key Milestones

| Milestone | Target Date | Criterion | Status |
|-----------|-------------|-----------|--------|
| **M0** Engine setup validated | Aug 15, 2026 | FreeCAD engine + scenario workflow ready | ✅ DONE |
| **M1** Phase 4 stable FSI | Dec 1, 2026 | Zero-order predictor validated, Aitken ω ∈ [0.03, 0.80] | 🔄 IN PROGRESS |
| **M2** Phase 5 complete (3 scenarios) | Jan 1, 2027 | WSS + pressure + displacement exported, plots done | 🔴 PLANNED |
| **M3** Phase 6 stabilization | Mar 15, 2027 | SUPG implemented, Neo-Hookean elastic wall works | 🔴 PLANNED |
| **M4** Phase 7 reports | Apr 15, 2027 | Interactive HTML + PDF dashboards ready | 🔴 PLANNED |
| **M5** 3D + MFEM ready | Jun 15, 2027 | 3D model imported, MFEM framework evaluated | 🔴 PLANNED |
| **M6** 🎯 DELIVERY | Sep 1, 2027 | Production FSI, full cardiac cycle, validation complete | 🔴 PLANNED |

---

## Technology Stack

| Layer | Tool | Version | Status | Notes |
|-------|------|---------|--------|-------|
| **Geometry / Parameterization** | FreeCAD | 1.0 | ✅ Active | Parametric aorta engine |
| **Mesh** | Gmsh | 4.13 | ✅ Active | 2D mesh generation, refinement |
| **CFD/FSI** | DolfinX | 0.11 | ✅ Active | Finite element framework (FEniCSx) |
| **Linear Algebra** | PETSc | 3.25 | ✅ Active | Solvers, preconditioners |
| **0D Models** | Python/NumPy/SciPy | 3.11 | ✅ Active | Windkessel, elastance models |
| **Data Processing** | Pandas / NumPy | current | ✅ Active | CSV export, analysis |
| **Visualization** | ParaView | 5.13 | ✅ Active | VTK rendering, videos |
| **CLI/Scripting** | Python | 3.10+ | ✅ Active | Automation, logging |
| **Production FSI** | MFEM | next | 🔴 Phase 8 | 3D monolithic solver |
| **Tissue Mechanics** | FEBio | next | 🔴 Phase 8 | Alternative to MFEM |
| **ML/Analytics** | TensorFlow/PyTorch | future | 🔴 Beyond scope | ROM, surrogate models |

**Environment:**
```
OS: Ubuntu 22.04 LTS
Dependencies: installed via conda/pip
Version control: Git + GitHub
CI/CD: GitHub Actions (optional)
HPC: MPI-enabled PETSc for Phase 8 scaling
```

---

## Risk Register

| Risk | Impact | Probability | Mitigation | Contingency |
|------|--------|-------------|------------|-------------|
| **FSI divergence** on pathological cases (Phase 5) | 🔴 High | Medium | Pre-test with simple stenosis; use smaller Δt (1e-4) | Revert to quasi-static analysis |
| **Aitken ω saturation** at 0.80 (slow convergence) | 🟡 Medium | Medium | Consider Anderson mixing as alternative | Hybrid Anderson-Aitken |
| **SUPG implementation** complexity (Phase 6) | 🟡 Medium | Low | Reference FEniCSx docs; test on 2D cylinder first | Use commercial SUPG library |
| **3D mesh quality** issues (Phase 8) | 🔴 High | Medium | Validate with industry CT scans early; use mesh metrics | Manual mesh refinement |
| **MFEM monolithic solver** instability | 🔴 High | Medium | Benchmark against DolfinX on 2D; staged migration | Fall back to FEBio |
| **Memory overflow** on 3D (100K DOF) | 🟡 Medium | Low | Use distributed memory (MPI); iterative solvers | Reduce mesh resolution |
| **Reproducibility issues** (code drift) | 🟡 Medium | Low | Pin all library versions; document env | Dockerfile for reproducibility |

---

## Recommended Actions (October 2026)

### Immediate (This Week)
- [ ] **Confirm Phase 4 FSI stability** on baseline case (400 steps)
- [ ] **Run 50-step test** with `clot_severe` BC (stress-test convergence)
- [ ] **Document convergence logs** (Aitken ω vs iteration count)
- [ ] **Verify mass balance** and energy conservation

### Next 2 Weeks (Mid-Oct)
- [ ] **Plan Phase 5 BC integration:**
  - [ ] Create adapter: `turbine_bc_*.npz` → inlet Q(t)
  - [ ] Test on 1 baseline case
  - [ ] Create multi-scenario harness (loop over 3 cases)
- [ ] **Set up CSV export pipeline:**
  - [ ] WSS per element (time-series)
  - [ ] Pressure at 5 monitoring points
  - [ ] Interface displacement (L2 norm)

### Next Month (Oct–Nov)
- [ ] **Complete Phase 5:**
  - [ ] Run baseline FSI (3 cycles)
  - [ ] Run hypertension FSI
  - [ ] Run stenosis FSI (with fallback Δt if needed)
  - [ ] Extract + plot comparisons
  - [ ] Draft clinical interpretation

### Strategic (Q1 2027)
- [ ] **Phase 6 planning:**
  - [ ] Allocate developer time for SUPG implementation
  - [ ] Review FEniCSx SUPG examples
  - [ ] Procure Neo-Hookean material model code
- [ ] **Phase 8 prep:**
  - [ ] Identify 3D aorta dataset (CT/synthetic)
  - [ ] Evaluate MFEM licensing + installation
  - [ ] Budget for HPC resources (CPU-hours)

---

## Project Metrics & KPIs

| Metric | Target | Current | Status |
|--------|--------|---------|--------|
| **Phase 4 stability** | < 10 divergences per 100 runs | TBD | 🔄 Testing |
| **FSI iterations/step** | 5–8 (Aitken) | TBD | 🔄 Profiling |
| **Solver wall-clock time** | < 3.0 s/step (2D) | TBD | 🔄 Benchmarking |
| **Mass conservation** | < 1% error | ~0.5% | ✅ Good |
| **Mesh independence** | < 2% error (V_peak) | TBD | 🔄 Phase 5 check |
| **Publication readiness** | Peer-review quality (Phase 7) | Not yet | 🔴 Planned |
| **3D scalability** | < 50 s/step (100K DOF) | N/A | 🔴 Phase 8 target |

---

## Code Organization (Key Directories)

```
Global-model-of-blood-circulation-and-pathology/
├── src/
│   ├── geometry/          # FreeCAD models, Gmsh scripts
│   │   └── build_schematic_2d.py
│   ├── solvers/           # CFD, FSI implementations
│   │   ├── ns_aorta_bdf2.py
│   │   ├── aorta_fsi_bdf2.py
│   │   └── windkessel_bc.py
│   ├── data/              # Clinical data parsing
│   │   ├── parser_who_gbd.py
│   │   └── export_turbine_bc.py
│   └── utils/             # Utilities (logging, I/O)
├── tests/                 # Unit + regression tests
├── results/               # Output: VTK, CSV, figures
├── docs/                  # Documentation, equations
└── ROADMAP_GANTT.md       # This file
```

---

## Glossary

| Term | Definition |
|------|-----------|
| **ALE** | Arbitrary Lagrangian-Eulerian; mesh follows deforming boundaries |
| **BDF2** | 2nd-order Backward Differentiation Formula (time integration) |
| **CFD** | Computational Fluid Dynamics |
| **DOF** | Degrees of Freedom (mesh nodes × variables) |
| **FSI** | Fluid-Structure Interaction |
| **IPCS** | Incremental Pressure Correction Scheme |
| **KDTree** | K-d spatial tree for nearest-neighbor queries |
| **MFEM** | Modular Finite Element Methods library |
| **SUPG** | Streamline Upwind Petrov-Galerkin (stabilization) |
| **WSS** | Wall Shear Stress (τ = μ ∂u/∂n at wall) |

---

## References & Resources

### Papers
- Quarteroni et al., *Computational Methods for FSI* (Springer)
- Takizawa & Tezduyar, *ALE methods for 3D FSI* (IJCFD)
- Peskin, *The Immersed Boundary Method* (Acta Numerica)

### Software Documentation
- **FEniCSx:** https://fenicsproject.org/
- **Gmsh:** https://gmsh.info/
- **PETSc:** https://petsc.org/
- **ParaView:** https://www.paraview.org/
- **MFEM:** https://mfem.org/

### Clinical Context
- WHO CVD guidelines (https://www.who.int/)
- GBD Disease Burden Study (https://www.healthdata.org/gbd)
- Arterial hemodynamics literature (Journals: AJP, Circulation)

---

**Status:** 🔄 On track for Q3 2027 delivery. Phase 4 (FSI) active; Phase 5 gates Phase 6 dependencies.

**Last Updated:** October 2026 | **Next Review:** November 1, 2026

