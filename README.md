# Global Model of Blood Circulation and Pathology

**Portfolio project:** Blood pump CFD and aortic FSI in FEniCSx. Numerical modeling of blood flow and fluid-structure interaction with FEniCSx/DolfinX.

---

## Overview

This portfolio project combines:
- **Computational Fluid Dynamics (CFD)** for blood pump analysis
- **Fluid-Structure Interaction (FSI)** for aortic wall modeling
- **0D-3D coupled systems** for realistic hemodynamics
- **Pathology scenarios** for disease modeling

---

## ✅ Completed

- **Steady-state 3D geometry analysis:** Comparison of blade designs (N1, N2, N3) and staggered rows with mesh-independence checks (pressure drop, hemolysis, torque)
- **Unsteady Navier-Stokes:** Validated on Kármán vortex street, extended to 3D pump
- **2D rigid-wall aorta model:**
  - 0D heart model at inlet
  - Windkessel outlets on all four branches
  - Clot simulation
  - Comparison of IPCS, Newton-CN, and Newton-BDF2 schemes
- **0D pathology scenarios:** Five cases (hypertension, diabetes, stenosis, combined) with boundary conditions exported for CFD

---

## 🔄 In Progress

- **Strongly coupled BDF2-ALE FSI:** Elastic wall, traction transfer, mesh velocity, GCL stability
- **"Balloon" test:** Smooth expansion of aortic arch under pulsatile pressure

---

## 📋 Planned

- Variable wall thickness and large-deformation wall models
- MFEM/FEBio coupling, 3D bifurcation and aortic arch with turbine
- Plaques, stents, rupture risk assessment, hemolysis analysis
- ML models for diagnostics and patient-state prediction

---

## Technologies

- **FEniCSx/DolfinX** – Finite element framework
- **Python** – Primary development language
- **CFD & FSI** – Cardiovascular flow modeling
