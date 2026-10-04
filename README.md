# Global Model of Blood Circulation and Pathology

**Portfolio project:** Blood pump CFD and aortic FSI in FEniCSx. Numerical modeling of blood flow and fluid-structure interaction with FEniCSx/DolfinX.

---

## 🎯 Project Status

**Active Timeline:** July 2026 → Q2 2027 | **Current Phase:** 🔄 Ph.4 FSI Solver

| Phase | Status | Description |
|-------|--------|-------------|
| Ph.0 | ✅ Done | Turbine Engine Development – 9 blade geometries (FreeCAD) |
| Ph.1 | ✅ Done | Parametric Geometry & Meshing – 2D aorta + 3 branches |
| Ph.2 | ✅ Done | Rigid CFD – BDF2 + Windkessel outlets, 3-scheme validation |
| Ph.3 | ✅ Done | Data Science – 20 NPZ clinical scenarios, heart engine, thrombosis DB |
| **Ph.4** | **🔄 Active** | **Strongly Coupled FSI (ALE BDF2)** – Elastic wall, GCL stability |
| Ph.5 | 🔴 Planned | Pathology Runs – Connect 20 NPZ scenarios, WSS comparison |
| Ph.6 | 🔴 Planned | SUPG Stabilization + 3D Neo-Hookean + 3D aortic arch |
| Ph.7 | 🔴 Planned | ML Health Prediction – Diagnostics & patient risk scoring |

---

## Overview

This portfolio project combines:
- **Computational Fluid Dynamics (CFD)** for blood pump analysis
- **Fluid-Structure Interaction (FSI)** for aortic wall modeling
- **0D-3D coupled systems** for realistic hemodynamics
- **Pathology scenarios** for disease modeling

---

## ✅ Completed Work

### Phase 0: Turbine Engine Development
- **Parametric generator:** FreeCAD Part API (no PartDesign) – fully scriptable
- **9-case blade sweep matrix:** N1, N2 (flat/stagger), N3 (flat/stagger) with coarse/medium/fine meshes
- **Mesh-independence study:** Torque, pressure drop, hemolysis across resolutions
- **Output:** 9 MSH v2 ASCII files ready for CFD

### Phase 1: Geometry & Meshing
- **2D rigid-wall aorta:** Parametric inlet, 3 branch outlets
- **Gmsh meshing pipeline:** FreeCAD STEP → Salome (surface repair) → Gmsh (physical groups) → MSH v2 ASCII
- **File:** `schematic_rects.msh`

### Phase 2: Rigid CFD
- **Unsteady Navier-Stokes:** Validated on Kármán vortex street, extended to 3D pump
- **Three time-stepping schemes:** IPCS, Newton-CN, Newton-BDF2 comparison
- **Boundary conditions:**
  - 0D heart model at inlet (time-varying elastance, valve diodes)
  - Windkessel RC outlets on all four branches
  - Wall clot simulation
- **Validation metrics:** Mass error < 1.5%, velocity leakage fix, skew-symmetric convection
- **Artificial stabilization:** Backflow control, adaptive time-stepping

### Phase 3: Data Science & Clinical Scenarios
- **Thrombosis database:** 12 nosologies (ESC 2023, AHA/ACC guidelines, WHO data)
- **Heart engine:** `heart_engine.py` with 15 thrombosis variants
- **20 NPZ boundary conditions:** Turbine inlet flows for 5 normal + 15 pathological cases
  - **Baseline:** Normal circulation
  - **Hypertension variants:** Mild, moderate, severe
  - **Metabolic:** Diabetes + vessel stiffness
  - **Stenosis:** Aortic valve, cholesterol-induced
  - **Combined:** Severe multifactor disease, "WHO killer combo"
  - **Critical events:** Heart attack, stroke, cardiac arrest
  - **Thrombosis:** Mild, moderate, severe clot, A-fib

---

## 🔄 In Progress: Phase 4

### Strongly Coupled BDF2-ALE FSI Solver

**Working components:**
- ✓ Submesh extraction (fluid/solid interface tracking)
- ✓ KDTree DOF mapping (~1e-9 m precision)
- ✓ Aitken underrelaxation (ω 0.03→0.6 adaptive)
- ✓ Zero-order predictor for stability
- ✓ VTK output every 10 FSI iterations (40 frames per run)
- ✓ Full convergence logging

**Current challenges:**
- ⚠ **Velocity oscillations** at high Reynolds number → SUPG stabilization needed
- ⚠ **Added-mass instability** (ρ_fluid ≈ ρ_solid = 1060 kg/m³) → mitigated with DT=2.5e-4
- ⚠ **Sharp corner singularities** in 2D aorta → will smooth in Ph.6
- ⚠ **Cantilever bending:** Side branches bend instead of radially expanding (2D limitation)

**Milestone gates:**
- 🎯 **Gate 4A:** 400 stable steps, no divergence → unlock Ph.5
- 🎯 **Gate 4B:** SUPG added, Vmax < 2 m/s, smooth velocity/pressure fields
- 🎯 **Gate 4C:** Wall displacement 0.3–1.0 mm, physiology confirmed

---

## 📊 Blade & Turbine Analysis

### 9-Case Mesh-Independence Study
Results show impact of blade count and arrangement on hemodynamics:
- **N1 (1 blade):** Baseline, smooth torque, low hemolysis
- **N2 flat (2 rows, parallel):** Doubled torque, moderate hemolysis
- **N2 stagger (chess pattern):** Chaotic flow, reduced torque, complex wake
- **N3 flat (3 rows, parallel):** Highest torque, maximum hemolysis risk – **mesh independence confirmed**
- **N3 stagger (3 rows, staggered):** Near-zero torque, flow redistribution – 239k fine-mesh exceeds 200k cell limit

### PrePoMax FEA (Rotor Mechanics)
- **Materials:** Aluminum rotor (E=70 GPa), steel casing (E=210 GPa)
- **Loading:** Pump pressure 13–40 kPa → peak elastic displacement ~27 microns
- **Centrifugal loading:** Via density card (high-speed rotor)

---

## 📋 Planned Features (Ph.5–7)

### Phase 5: Pathological Runs
- Connect 20 NPZ BCs → FSI inlet, run all scenarios
- Compare WSS, pressure, displacement across healthy vs. pathology
- Export hemodynamic features for ML training

### Phase 6: Advanced Modeling
- **SUPG stabilization:** Eliminate velocity noise
- **Neo-Hookean wall:** Large-deformation FSI (radial expansion vs. bending)
- **3D aortic geometry:** Smooth arch, bifurcation, realistic branch angles
- **Turbulence modeling:** k-ε or k-ω SST for realistic Re

### Phase 7: ML & Diagnostics
- Train models on hemodynamic features (WSS, pressure, displacement)
- Predict patient risk from CFD outputs
- Image recognition for diagnostic assistance

---

## Technologies

| Component | Technology | Status |
|-----------|-----------|--------|
| **Parametric CAD** | FreeCAD 1.0 (Part API) | ✅ Active |
| **Geometry repair** | Salome 9.x | ✅ Active |
| **Meshing** | Gmsh 4.x | ✅ Active |
| **CFD/FSI solver** | FEniCSx 0.11 / DolfinX | ✅ Active |
| **Linear algebra** | PETSc 3.25 + SLEPc | ✅ Active |
| **Scientific computing** | NumPy, SciPy, Matplotlib | ✅ Active |
| **Visualization** | ParaView + PyVista | ✅ Active |
| **Advanced FSI** | MFEM / FEBio | 🔴 Planned |
| **Machine learning** | PyTorch / TensorFlow | 🔴 Planned |

---

## Risk Register

| Risk | Severity | Mitigation |
|------|----------|-----------|
| Velocity oscillations (high Re) | 🔴 Critical | Add SUPG stabilization; Petrov–Galerkin for convection |
| Added-mass instability | 🔴 Critical | Zero-order predictor; DT=2.5e-4; monitor >200 steps |
| N3 stagger fine mesh (239k cells) | ⚠ Warning | Reduce Gmsh size factor; stay <200k limit |
| 2D cantilever bending | ⚠ Warning | Accept for Ph.4; smooth 3D geometry in Ph.6 |
| Ubuntu 26.04 env setup | ⚠ Warning | `pip install meshio pandas`; register fenicsx kernel |

---

## Next Actions

### 🔥 Immediate (This Week)
1. **Fix Ubuntu environment:** `pip install meshio pandas` + kernel registration
2. **Add SUPG to aorta_fsi_bdf2.py:** Eliminate velocity noise
3. **Verify 400-step FSI run:** Check displacement 0.3–1.0 mm, Vmax < 2 m/s

### 📅 Short-term (October 2026)
4. Connect `turbine_bc_0_baseline.npz` as FSI inlet BC (replace hardcoded velocity)
5. Run `turbine_bc_7_clot_severe.npz` scenario, export WSS comparison
6. Fix N3-stagger-fine mesh (239k → <200k cells)

### 🗓 Medium-term (Nov–Dec 2026)
7. Unsteady turbine solver with Reynolds turbulence model
8. Neo-Hookean wall model + 3D aorta (radial expansion validation)
9. Kármán vortex street: complete STEP→Salome→Gmsh→CFD pipeline

---

## File Structure

```
├── turbine_generator/          # Ph.0: Blade generator
│   ├── config.py               # TurbineConfig dataclass
│   ├── blade.py, hub.py, ...   # Geometry components
│   ├── generate.py             # Main entry point
│   └── README.md               # Full API docs
├── 0D_models/                  # Ph.3: Clinical scenarios
│   ├── heart_engine.py         # Elastance + Windkessel
│   ├── heart_engine_sweep.py   # 20-scenario automation
│   └── thrombosis_*.py         # Disease variants
├── 0D_results/                 # Output: 20 NPZ files
│   ├── turbine_bc_normal.npz
│   ├── turbine_bc_0_baseline.npz
│   ├── turbine_bc_1_hypertension*.npz
│   ├── turbine_bc_3_aortic_stenosis.npz
│   ├── turbine_bc_5_heart_attack.npz
│   └── turbine_bc_clot_*.npz
├── CFD_rigid/                  # Ph.2: Rigid solver
│   ├── ns_staggered_plates_pulsatile_newton.py
│   ├── run_meshindep_sweep.ipynb
│   └── build_plate2d_*.py      # Scheme documentation
├── CFD_FSI/                    # Ph.4: FSI solver (active)
│   ├── aorta_fsi_bdf2.py       # Main FSI loop
│   ├── fsi_runner.ipynb
│   └── validate_interface.py
├── meshes/                     # Gmsh outputs
│   ├── schematic_rects.msh     # 2D aorta
│   ├── plate2d.msh
│   └── turbine_*.msh           # 9 blade cases
├── project_dashboard_v2.html   # Interactive dashboard
└── README.md                   # This file
```

---

## Quick Links

- **Dashboard:** `project_dashboard_v2.html` (open in browser for interactive Gantt chart, risk register, and next actions)
- **0D scenarios:** `0D_results/turbine_bc_*.npz` (NumPy binary format: t, Q_AO, u_mean, P_AO, T_PERIOD)
- **Mesh pipeline:** FreeCAD → Salome → Gmsh → DolfinX
- **Visualization:** ParaView for VTK, Matplotlib for publication figures

---

## References

- FEniCSx Documentation: https://fenicsproject.org/
- ESC 2023 Guidelines: European Society of Cardiology
- PETSc/SLEPc: https://petsc.org/
- Gmsh: https://gmsh.info/

---

**Last Updated:** October 2026  
**Contact:** mircool6  
**License:** Portfolio / Educational Use
