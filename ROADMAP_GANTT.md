[masterplan.txt](https://github.com/user-attachments/files/33024169/masterplan.txt)
================================================================================
MASTERPLAN — Blood Pump CFD & Aortic FSI Simulation
================================================================================
Project Start : August 2026
Last Updated  : October 2026
Stack         : Python · DolfinX 0.11 · PETSc · Gmsh · ParaView
================================================================================


--------------------------------------------------------------------------------
PHASE 1 — Geometry & Mesh  [DONE — August 2026]
--------------------------------------------------------------------------------
Goal: Build a 2D schematic aorta geometry with branches and label boundaries.

Files:
  build_schematic_2d.py         Gmsh: aorta + 3 branches, rectangular scheme
  schematic_2d.msh              Base mesh (rigid walls)
  schematic_rects.py            Refined rectangular geometry
  schematic_rects.msh           Mesh for rigid CFD

Result: Parametric 2D mesh with inlet, outlets and wall tags.
        Preview saved to preview.png


--------------------------------------------------------------------------------
PHASE 2 — Rigid CFD  [DONE — August–September 2026]
--------------------------------------------------------------------------------
Goal: Solve unsteady Navier-Stokes on a rigid aorta, validate solvers.

Files:
  ns_aorta_4outlet_heart_pulsatile.py         IPCS pulsatile flow, 4 outlets
  ns_aorta_4outlet_heart_pulsatile_RC_fixed.py RC outlets (Windkessel)
  compare_aorta_3schemes.py                   IPCS / Newton-CN / Newton-BDF2
  compare_aorta_ipcs_newton.py                Detailed two-scheme comparison
  heart_engine.py                             0D heart model (time-varying elastance)
  aorta_bdf2.py                               BDF2 solver (final scheme)

Key IPCS fixes applied:
  - Velocity leakage fix (apply_lifting)
  - Skew-symmetric convection
  - Artificial viscosity + backflow stabilization
  - Adaptive dt + mass error < 1.5%
  - Windkessel synchronization


--------------------------------------------------------------------------------
PHASE 3 — 0D Pathophysiology & Boundary Conditions  [DONE — September 2026]
--------------------------------------------------------------------------------
Goal: Generate physiologically grounded BCs for 12 clinical scenarios.

Files:
  heart_engine_clot.py          0D model with thrombosis
  heart_engine_sweep.py         Parameter sweep across 5 base scenarios
  thrombosis_data_fetcher.py    ESC/AHA/WHO data (12 nosologies)
  thrombosis_summary.md         Top-5 most lethal conditions
  markov_cfd_analysis.py        Markov-chain progression analysis

Exported NPZ files (0D_results/):
  turbine_bc_normal.npz               Normal / baseline
  turbine_bc_0_baseline.npz           Baseline (labeled)
  turbine_bc_1_hypertension.npz       Hypertension
  turbine_bc_1_hypertension_mild.npz  Mild hypertension
  turbine_bc_2_diabetes_stiff.npz     Diabetes (arterial stiffness)
  turbine_bc_3_aortic_stenosis.npz    Aortic stenosis
  turbine_bc_3_cholesterol_stenosis.npz  Cholesterol-driven stenosis
  turbine_bc_4_severe_combined.npz    Severe combined pathology
  turbine_bc_4_who_killer_combo.npz   WHO killer combination
  turbine_bc_5_heart_attack.npz       Myocardial infarction
  turbine_bc_6_stroke.npz             Ischemic stroke
  turbine_bc_7_cardiac_arrest.npz     Cardiac arrest
  turbine_bc_clot_mild.npz            Mild thrombus
  turbine_bc_clot_moderate.npz        Moderate thrombus
  turbine_bc_clot_severe.npz          Severe thrombus
  turbine_bc_clot_afib.npz            Atrial fibrillation + thrombus

PDF dashboards:
  dashboard_0_baseline.pdf
  dashboard_4_who_killer_combo.pdf
  0D_biophysical_sweep.pdf


--------------------------------------------------------------------------------
PHASE 4 — FSI (Fluid-Structure Interaction)  [IN PROGRESS — Sep–Oct 2026]
--------------------------------------------------------------------------------
Goal: Strongly-coupled partitioned ALE FSI on a 2D aorta with elastic walls.

Solver architecture:
  Fluid  (BDF2 ALE Navier-Stokes)
      |  traction  sigma*n
  Interface DOF mapping  (KDTree, error < 1e-9 m)
      |  displacement  d
  Solid  (2D Linear Elasticity, Plane Strain)
      |
  ALE mesh smoothing  (Laplace)
      |
  Aitken relaxation  (dynamic omega)

Files:
  schematic_rects_fsi.py        FSI mesh (fluid=100, wall=200, iface=6)
  schematic_rects_fsi.msh       3683 nodes, 7792 elements
  aorta_fsi_bdf2.py             Main FSI solver
  FSI_DEVELOPMENT_STATUS.md     Detailed problem/fix log
  FSI_Scheme_Description.md     Mathematical description of the scheme

What works:
  [x] Fluid/solid submesh extraction from a single mesh
  [x] Interface DOF mapping via KDTree (error < 1e-9 m)
  [x] FSI loop convergence (Aitken, omega from 0.03 up to 0.6+)
  [x] Zero-order predictor (eliminates artificial velocity shocks)
  [x] VTK export every 10 steps  ->  40 frames per 400 steps
  [x] Detailed log: step / fsi-iter / traction / pressure / displacement

Current parameters:
  DT                 = 2.5e-4 s
  MAX_STEPS          = 400
  WRITE_EVERY        = 10
  FSI_RELAX_MAX      = 0.80
  MAX_INTERFACE_DISP = 2.0e-3  (2 mm)

Known limitations:
  - 2D Plane Strain: side branches bend like cantilevers instead of radially
    expanding. Accepted as working model for this phase.
  - Numerical noise in velocity field at high Re -> needs SUPG stabilization.
  - Sharp branch corners cause velocity singularities.


--------------------------------------------------------------------------------
PHASE 5 — Pathological Runs & Analysis  [PLANNED]
--------------------------------------------------------------------------------
  [ ] Connect turbine_bc_*.npz files as inlet BCs in aorta_fsi_bdf2.py
  [ ] Run: baseline, hypertension, clot_severe (minimum 3 scenarios)
  [ ] Export WSS, pressure, wall displacement to CSV
  [ ] Build comparison plots (matplotlib / PDF report)
  [ ] Write clinical interpretation for each scenario


--------------------------------------------------------------------------------
PHASE 6 — Stabilization & Quality Improvement  [PLANNED]
--------------------------------------------------------------------------------
  [ ] SUPG stabilization for the convective term in Navier-Stokes
  [ ] Nonlinear wall mechanics (Neo-Hookean instead of linear elasticity)
  [ ] Smoother branch geometry (rounded corners, no singularities)
  [ ] 3D geometry: upgrade from 2D cross-section to full 3D aorta
  [ ] Migration to MFEM / FEBio for production FSI (monolithic solver)


--------------------------------------------------------------------------------
PHASE 7 — Visualization & Reporting  [PLANNED]
--------------------------------------------------------------------------------
  [ ] Interactive dashboard (generative UI / HTML5)
  [ ] PDF report: scheme description + results per scenario
  [ ] ParaView animations for each clinical scenario
  [ ] Comparison table: healthy aorta vs pathologies


================================================================================
TECHNOLOGY STACK
================================================================================
  Mesh generation   Gmsh                    ACTIVE
  CFD solver        DolfinX 0.11 / FEniCSx  ACTIVE
  Linear algebra    PETSc 3.25              ACTIVE
  0D models         Python / NumPy / SciPy  ACTIVE
  Visualization     ParaView + Matplotlib   ACTIVE
  FSI (production)  MFEM / FEBio            NOT STARTED
  Machine learning  TensorFlow / PyTorch    NOT STARTED
================================================================================
