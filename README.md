# Global-model-of-blood-circulation-and-pathology
Portfolio project: blood pump CFD and aortic FSI in FEniCSx. Done: 3D blade comparison, unsteady Navier-Stokes, 2D rigid aorta with 0D heart and Windkessel outlets, pathology scenarios. In progress: strongly coupled BDF2-ALE FSI. Planned: plaques, stents, hemolysis, ML diagnostics.

Blood Pump CFD & Healthcare Simulation

Portfolio project for Siemens Healthineers: numerical modeling of blood flow and FSI with FEniCSx/DolfinX.

Completed
Steady-state analysis of 3D geometries: comparison of blade designs N1/N2/N3 and staggered rows, mesh-independence checks (dP, hemolysis, torque).
Unsteady Navier-Stokes: validated on a Kármán vortex street, then extended to the 3D pump.
2D rigid-wall aorta: 0D heart model at the inlet, Windkessel outlets on all four branches, a clot, and a comparison of IPCS, Newton-CN and Newton-BDF2.
0D pathology scenarios: five cases (hypertension, diabetes, stenosis, combined) with boundary conditions exported for CFD.
In Progress
Strongly coupled BDF2-ALE FSI: elastic wall, traction transfer, mesh velocity, GCL.
"Balloon" test: smooth expansion of the aortic arch under pulsatile pressure.
Planned
Variable wall thickness and large-deformation wall model.
MFEM/FEBio coupling, 3D bifurcation and aortic arch with a turbine.
Plaques, stent, rupture risk, hemolysis analysis.
ML models for diagnostics and patient-state prediction.
