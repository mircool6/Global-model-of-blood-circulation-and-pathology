MASTER PLAN - BLOOD PUMP CFD & HEALTHCARE SIMULATION (Siemens Portfolio Project)
MAIN AND COMPLETED PHASES (1-3)

1. Mesher improvement

[DONE] Switched to direct generation via gmsh_from_step.py instead of salome_mesher.py.

2. Steady-state geometry analysis (Phases 1 and 2)

[DONE] Comparison of blades N1, N2, N3 and stagger (1, 2, 3 blade rows).
[DONE] Testing of fine, middle and coarse meshes (limit 200k cells).
[DONE] 9 plots for N1, N2, N3, 6 plots for N3/stagger, 3 hemodynamics tables.

3. Unsteady analysis and Kármán vortex street (Phase 3)

[DONE] Kármán vortex street implemented (FreeCAD script -> Gmsh).
[DONE] Moved from steady-state to unsteady Navier-Stokes (implicit Backward Euler), animation run.
CFD solvers used and their status
Method	Dimension / case	Time and nonlinear scheme	Status and conclusion
Stationary Stokes -> stationary Newton	3D blade/turbine geometry sweep, Phases 1-2	Linear mixed Stokes as initial guess (direct LU/MUMPS) -> stationary mixed Navier-Stokes, PETSc SNES Newton (newtonls, basic line search, direct LU/MUMPS)	[DONE] Used for N1/N2/N3, stagger and mesh-independence comparisons; obtained dP, hemolysis, stagnation and torque.
Picard (frozen convection)	2D benchmark: cylinder / Kármán vortex street	Backward Euler, usually 1-2 Picard iterations per step	[DONE] Validation of the 2D solver and Kármán vortices obtained.
Picard (frozen convection)	3D: turbine / blade transient cases	Backward Euler, 2 Picard iterations, grad-div stabilization	[DONE] Working early unsteady solver; used for mesh, CFL and pressure-drop tests.
IPCS / pressure-correction (fixed)	2D rigid aorta, 4 outlets + RC/Windkessel	Backward Euler: tentative velocity -> Poisson pressure correction -> velocity correction	[DONE - STABLE] Critical bugs fixed (see Phase 4). mass_rel < 1.5%, du_proj < 0.005.
Monolithic mixed Newton	2D rigid aorta, 4 outlets + RCR/Windkessel	Backward Euler, coupled P2/P1 solution of velocity and pressure; PETSc SNES Newton (newtonls)	[DONE] Used as the reference for comparison with IPCS.
HEMODYNAMICS AND AORTA (PHASE 4) - COMPLETED

4. Blood flow modeling in the aorta (rigid-wall CFD)

[DONE] Base aorta geometry (inlet -> ascending -> arch -> descending, branches).
[DONE] DolfinX set up for aortic flow simulation.
[DONE] 0D pulsating heart model connected at the inlet.
[DONE] RC/RCR Windkessel on all 4 outlets.
[DONE] Clot: static blocking region in the aorta, USE_CLOT=True, Re=500.
[DONE] Comparison of 3 schemes: IPCS vs Newton-CN vs Newton-BDF2 (aorta_compare.py).
[DONE] Pacemaker (0D Heart Engine) - heart_engine_clot.py.

Critical IPCS fixes (all done)

[DONE] Fix A3: assemble_matrix(A3, a3, bcs=bcu3) - removed velocity leakage through the walls (main bug).
[DONE] Fix apply_lifting in Step 3 for correct Dirichlet BCs.
[DONE] Skew-symmetric convection in F1 - suppressed energy generation when div(u) != 0.
[DONE] Artificial viscosity (mu_art = 0.1 * U_IN * hmin).
[DONE] Backflow stabilization at the outlets.
[DONE] Adaptive dt with ABORT (V_ABORT = 10 * vmax_prev, not hardcoded).
[DONE] Windkessel synchronized with the actual t_now += dt_try.
[DONE] Removed the manual set_bc(uh, bcu) after projection - mass error now < 1.5%.
[DONE] Detailed diagnostics: du_star, du_proj, max|phi|, Qin*, dQin, flux_defect.
PATHOPHYSIOLOGY AND 0D MODELS (PHASE 5) - CURRENT STAGE

5. Knowledge base on thrombosis and risk factors

[DONE] thrombosis_stats.csv - 12 disease entities (baseline parameters for calibration).
[DONE] WHO GBD / UCI / PubMed data parser (automatic collection of CVD statistics and risk factors).
[DONE] cfd_biophysical_triggers.csv - biophysical thresholds for CFD (Low WSS < 0.4 Pa, High WSS > 50 Pa, sugar, cholesterol).
[DONE] who_cvd_risk_factors.csv - global mortality by risk factor (what destroys the artery).
[DONE] Data visualization (thrombosis_report.pdf) - automatic plot generation via matplotlib.

6. 0D pathological sweep (heart_engine_sweep.py)

[DONE] Mathematical core of the curve simulation (engine physics):
Heart: Time-Varying Elastance E(t) (ODE integration: P_LV = E(t)*(V_LV - V0)).
Vessels: 0D Windkessel (RC chains, where Compliance C = elasticity, Resistance R = periphery).
Valves: diode logic based on the pressure gradient (opens when P_LV > P_AO).
[DONE] heart_engine_sweep.py - 5 scenarios based on WHO biophysics:
0_baseline: healthy patient (normal C and R).
1_hypertension: hypertension (high systemic resistance R_sys, High WSS trigger).
2_diabetes_stiff: diabetes (reduced aortic compliance C_ao, glycocalyx damage, very large pulse pressure).
3_cholesterol_stenosis: stenosis/plaque (critical increase of aortic root resistance R_aortic).
4_who_killer_combo: WHO "killer" combination (hypertension + diabetes + plaque at once).
[DONE] Export of turbine_bc_<scenario>.npz for 2D/3D CFD (t, Q_AO, u_mean, P_AO, V_LV, P_LV).
[DONE] Visualization (0D_biophysical_sweep.pdf): Q(t), P(t), velocity profiles u(r) and combined ventricular PV loops.

7. 0D -> 2D CFD synchronization

[NOT DONE] Use Q_AO(t) and P_AO(t) from the .npz files as BCs for the IPCS solver.
[NOT DONE] Compare velocity and wall shear stress (WSS) fields in the 2D aorta for all 5 WHO scenarios.
[DONE] Inlet velocity profile: power law (1/7 Power Law) implemented for turbulent/transitional flow.
[NOT DONE] Exact Womersley inlet profile (if analytics are needed instead of 1/7).
[NOT DONE] Pulsatility index PI = (Q_max - Q_min) / Q_mean.
FSI (FLUID-STRUCTURE INTERACTION) (PHASE 6)
Resource budget and performance strategy
[DONE] Smoke-test budget adopted for the current DolfinX/FEniCSx prototype: DT = 2.5e-5 s, MAX_STEPS = 400, physical window = 0.010 s. Target: about 15 minutes of CPU time. This checks BDF2 convergence, CFL, mass balance, ALE/FSI coupling and field output; it is not a full cardiac cycle.
[RULE] Do not run the full period T = 0.8333 s (33,333 steps at this DT) until profiling is done and time/energy budget is agreed separately.
[PLAN] Optimization in order:
Profile DolfinX/PETSc: assembly, LU, number of Newton/FSI iterations, I/O and checkpoint frequency.
Speed up the algorithm without changing language: better preconditioning, MPI, adaptive output, shortened tests and verification on a short window.
For production/HPC, move to C++ (MFEM) with PETSc/Hypre; hot kernels can use C/C++ and a GPU path.
Fortran is acceptable for standalone numerical kernels/0D models, but is not the first path for MFEM/FEBio integration.
Use JavaScript only for UI, reports and interactive visualization; not for the CFD/FSI solver core.

8. Immersed FSI (combining MFEM and FEBio)

[REFERENCE] MLG_final_.pdf (University Koblenz, 2024): staged FSI/shape-optimization route.
[NOT DONE] Set up FEBio for simulating hyperelastic soft tissue (vessel wall).
[NOT DONE] Set up MFEM / MFEMiFSI for monolithic FSI coupling.
[NOT DONE] Benchmark 1: 2D pulsating channel with an elastic wall.
[NOT DONE] Benchmark 2: 3D carotid bifurcation with deformable walls.
[NOT DONE] Benchmark 3: 3D aortic arch FSI with turbine and elastic walls.
[NOT DONE] Add the turbine (blood pump) inside the aorta.
[NOT DONE] Unsteady pulsatile-flow run in the aorta with a working turbine.
[DONE] Separate 2D conforming FSI mesh created: schematic_rects_fsi.msh: fluid=100, wall=200, fsi_interface=6.
[IN PROGRESS] New strongly coupled BDF2-ALE FSI solver: traction transfer, elastic wall, mesh velocity, GCL and FSI residual convergence.
Current FSI task: physiological swelling of the aortic wall
[IN PROGRESS] FSI enabled on the whole fluid-wall interface, not only on a separate segment or in the bifurcation zone.
[NOT DONE] Geometrically nonlinear wall model (large deformation), so that pressure and traction themselves change the shape of the wall and ALE mesh.
[NOT DONE] Spatially varying wall thickness h(s), including the aortic arch and transitions to branches; thickness must enter the material stiffness.
[NOT DONE] Mechanics boundary conditions: clamp only the end sections (inlet/outlets), do not fix the whole outer wall surface.
[NOT DONE] "Balloon" smoke test: smooth normal expansion of the arch under pulsatile pressure, without cantilever bending of the branches as the main response.
[NOT DONE] Balloon-test diagnostics: max(d·n), change of local diameter/area, pressure, FSI residual, GCL and mass balance.
COMPLEX PATHOLOGIES, SYSTEMS AND ML (PHASES 7-9)

9. Global circulation model and pathologies

[NOT DONE] Create a global circulation scheme (coupling 3D and 0D models of the human body).
[NOT DONE] Modeling of cholesterol plaques and occlusion.
[NOT DONE] Simulation of aortic rupture risk.
[NOT DONE] Stent placement modeling.
[NOT DONE] Hemolysis analysis (calculate_metrics.py, Bludszuweit/Giersiepen-style blood damage model).

10. Machine learning in healthcare

[DONE] Dataset for the ML part downloaded from Kaggle.
[NOT DONE] ML model for predicting patient health (Siemens / Apple Watch approach).
[NOT DONE] ML assistant to help doctors with diagnosis (image recognition).
