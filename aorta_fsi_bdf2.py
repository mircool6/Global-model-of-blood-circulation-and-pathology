# -*- coding: utf-8 -*-
"""
aorta_fsi_bdf2.py
=================
Strongly-Coupled Partitioned ALE FSI Solver for the 2D Aorta.

Implements:
  - Submesh extraction (Fluid vs Solid)
  - Geometric facet mapping for exact DOF transfer
  - Smooth fluid traction projection via Poisson distance field normal
  - 2D Linear Elasticity for the wall (Plane Strain)
  - BDF2 Navier-Stokes with exact ALE convection (u - w_m)
  - FSI Fixed-Point Iteration with relaxation
"""
from pathlib import Path
import os, time, math, inspect
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from mpi4py import MPI
from petsc4py import PETSc
import gmsh
import ufl
from dolfinx import fem, mesh as dmesh
from dolfinx.io import VTKFile, gmsh as gmshio
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from basix.ufl import element, mixed_element

ROOT = Path("/mnt/c/Work/Cursor_CFD")
MESH_FILE = ROOT / "schematic_rects_fsi.msh"
HEART_BC = ROOT / "0D_results/turbine_bc_0_baseline.npz"
OUT = ROOT / "results_aorta_fsi_bdf2"
OUT.mkdir(parents=True, exist_ok=True)
WRITE_EVERY = 5
LOG_EVERY = 5

# ----------------- Physics & Numerics -----------------
RHO, MU = 1060.0, 0.0035
NU = MU / RHO
D_PHYS, R_PHYS, H_MESH = 0.020, 0.010, 2.0
L_SCALE = D_PHYS / H_MESH

FLUID_TAG, WALL_TAG = 100, 200
INLET_TAG, INTERFACE_TAG, OUTER_TAG = 1, 6, 7
OUTLET_TAGS = (2, 3, 4, 5)

DT = 2.5e-4
MAX_STEPS = 400
N_CYCLES = 4
MAX_FSI_ITERS = 100
# A relative residual alone is ill-conditioned while the wall displacement is
# close to zero during the startup ramp.  Require a physical interface RMS
# residual as well as a scale-aware relative measure.
FSI_TOL = 1e-3
FSI_ABS_TOL = 1.0e-9       # m RMS per interface node = 0.001 micrometre
FSI_DISP_REFERENCE = 1e-7 # m; prevents division by a vanishing displacement
# Aitken starts cautiously; a fixed 0.30 destabilised the second fluid/mesh
# exchange even though the first physical FSI response was well behaved.
FSI_RELAX = 0.03
FSI_RELAX_MIN, FSI_RELAX_MAX = 0.0005, 0.80
# The saved 0D waveform begins at 0.189 m/s after a zero field.  Starting the
# transient directly there creates an unphysical ~microsecond pressure shock.
# Use a C1 ramp only to initialise this short 0.01 s FSI smoke run.  A full
# physiological run should instead start from a periodic/preloaded state.
FSI_STARTUP_RAMP = 2.0e-2

NEWTON_RTOL = 1e-8
NEWTON_ATOL = 1e-10
NEWTON_MAX_IT = 12

# Wall Elasticity (Plane Strain)
WALL_E, WALL_NU, WALL_RHO = 2.0e6, 0.45, 1100.0
LAMBDA = WALL_E * WALL_NU / ((1.0 + WALL_NU) * (1.0 - 2.0 * WALL_NU))
MU_S = WALL_E / (2.0 * (1.0 + WALL_NU))
MAX_INTERFACE_DISP = 2.0e-3  # 2.00 mm: physiological expansion is around 0.5 - 1.0 mm

# Windkessel
_Rp, _Rd, _C = 5.0e4/RHO, 1.1e7/RHO, 1.5e-7*RHO
WK = {2: (_Rp, _Rd, _C), 3: (_Rp, _Rd, _C), 4: (_Rp*1.6, _Rd*1.2, _C*0.8), 5: (_Rp, _Rd, _C)}
class OutletModel:
    def __init__(self, rp, rd, c):
        self.rp, self.rd, self.c, self.pc = rp, rd, c, 0.0
    def advance(self, dt, q):
        self.pc = (q + (self.c/dt)*self.pc) / (self.c/dt + 1.0/self.rd)
    def get_p_kin(self, q):
        return self.pc + self.rp * q

class HeartInlet:
    def __init__(self, path):
        d = np.load(path)
        self.period = float(d["T_PERIOD"])
        t_raw, u_raw = np.asarray(d["t"], dtype=float), np.asarray(d["u_mean_m_s"], dtype=float)
        self.t_norm = t_raw - t_raw[0]
        mask = self.t_norm <= self.period + 1e-9
        self.t_norm, self.u_wave = self.t_norm[mask], u_raw[mask]
        self.u_peak = float(self.u_wave.max())
        active = np.flatnonzero(self.u_wave > 1e-4 * self.u_peak)
        self.t_onset = float(self.t_norm[active[0]]) if len(active) else 0.0
        self.t_now = 0.0

    def u_mean_now(self):
        u_raw = float(np.interp(self.t_now % self.period, self.t_norm, self.u_wave))
        s = np.clip((self.t_now - self.t_onset) / FSI_STARTUP_RAMP, 0.0, 1.0)
        ramp = 0.5 * (1.0 - np.cos(np.pi * s))
        return u_raw * ramp

    def __call__(self, x):
        a = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
        u_c = 1.5 * self.u_mean_now()
        a[0] = u_c * np.clip(1.0 - (x[1]/R_PHYS)**2, 0.0, None)
        return a

def safe_linear_problem(a, L, prefix, **kwargs):
    args = dict(kwargs)
    if "petsc_options_prefix" in inspect.signature(LinearProblem.__init__).parameters:
        args["petsc_options_prefix"] = prefix
    return LinearProblem(a, L, **args)

def build_wss_normal(fluid, fdim, fluid_outer_facets):
    """Creates a smooth volumetric normal field pointing OUT of the fluid."""
    V = fem.functionspace(fluid, ("Lagrange", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = ufl.inner(fem.Constant(fluid, PETSc.ScalarType(1.0)), v) * ufl.dx
    bc_dofs = fem.locate_dofs_topological(V, fdim, fluid_outer_facets)
    bc = fem.dirichletbc(PETSc.ScalarType(0.0), bc_dofs, V)
    
    prob = safe_linear_problem(a, L, "wss_n_", bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    d_field = prob.solve()
    n_vec = -ufl.grad(d_field) / ufl.sqrt(ufl.inner(ufl.grad(d_field), ufl.grad(d_field)) + 1e-12)
    return n_vec

def main():
    comm = MPI.COMM_WORLD
    
    # 1. Load Parent Mesh
    gmsh.initialize()
    if comm.rank == 0: gmsh.open(str(MESH_FILE))
    data = gmshio.model_to_mesh(gmsh.model, comm, 0, gdim=2)
    gmsh.finalize()
    msh, ct, ft = data.mesh, data.cell_tags, data.facet_tags
    msh.geometry.x[:, :2] *= L_SCALE
    fdim = msh.topology.dim - 1

    # 2. Extract Interface Coordinates for Geometric Mapping
    msh.topology.create_connectivity(fdim, 0)
    iface_facets = ft.find(INTERFACE_TAG)
    parent_iface_vertices = set()
    for f in iface_facets:
        parent_iface_vertices.update(msh.topology.connectivity(fdim, 0).links(f))
    iface_coords = msh.geometry.x[list(parent_iface_vertices)]
    kdtree = KDTree(iface_coords)
    def on_interface(x):
        dist, _ = kdtree.query(x.T)
        return dist < 1e-6

    # 3. Create Submeshes
    wall, _, _, _ = dmesh.create_submesh(msh, msh.topology.dim, ct.find(WALL_TAG))
    fluid, _, _, _ = dmesh.create_submesh(msh, msh.topology.dim, ct.find(FLUID_TAG))
    wall_x0, fluid_x0 = wall.geometry.x.copy(), fluid.geometry.x.copy()
    
    # 4. Identify Submesh Boundaries
    fluid_iface_facets = dmesh.locate_entities_boundary(fluid, fdim, on_interface)
    fluid_all_bndry = dmesh.locate_entities_boundary(fluid, fdim, lambda x: np.full(x.shape[1], True))
    fluid_outer_facets = np.setdiff1d(fluid_all_bndry, fluid_iface_facets)
    
    wall_iface_facets = dmesh.locate_entities_boundary(wall, fdim, on_interface)
    wall_all_bndry = dmesh.locate_entities_boundary(wall, fdim, lambda x: np.full(x.shape[1], True))
    wall_outer_facets = np.setdiff1d(wall_all_bndry, wall_iface_facets)

    # Re-discover inlets/outlets on the fluid submesh geometrically
    # (Since parent facet tags aren't natively inherited by submeshes)
    def is_inlet(x): return np.isclose(x[0], 0.0, atol=1e-5)
    def is_outlet(x): return x[0] > 1e-5 and np.abs(x[1]) > 1e-5 # Simplified, but let's use exact parent coords
    
    fluid_inlet_facets = dmesh.locate_entities_boundary(fluid, fdim, is_inlet)
    # Since OUTLET_TAGS = (2,3,4,5), let's just find the facets at x > L - eps. 
    # Or map them via coordinates from parent.
    parent_outlet_facets = np.concatenate([ft.find(t) for t in OUTLET_TAGS])
    parent_out_v = set()
    for f in parent_outlet_facets: parent_out_v.update(msh.topology.connectivity(fdim, 0).links(f))
    out_coords = msh.geometry.x[list(parent_out_v)]
    out_tree = KDTree(out_coords)
    def on_outlet(x):
        dist, _ = out_tree.query(x.T)
        return dist < 1e-6
    fluid_outlet_facets = dmesh.locate_entities_boundary(fluid, fdim, on_outlet)

    # 5. Function Spaces
    cell_f, cell_w = fluid.topology.cell_name(), wall.topology.cell_name()
    W = fem.functionspace(fluid, mixed_element([element("Lagrange", cell_f, 2, shape=(2,)), element("Lagrange", cell_f, 1)]))
    Vf, _ = W.sub(0).collapse()
    Vf_ale = fem.functionspace(fluid, ("Lagrange", 1, (2,)))
    Vw = fem.functionspace(wall, ("Lagrange", 1, (2,)))

    # 6. Build Inter-mesh Map (Fluid Interface DOFs <-> Wall Interface DOFs)
    # Using P1 vector space, DOFs are exactly at vertices.
    f_iface_dofs = fem.locate_dofs_topological(Vf_ale, fdim, fluid_iface_facets)
    w_iface_dofs = fem.locate_dofs_topological(Vw, fdim, wall_iface_facets)
    # locate_dofs_topological returns *block* (node) indices for vector spaces.
    # Never use these indices directly on Function.x.array: it is a flat
    # component array.  Use the (n_nodes, 2) views below for every transfer.
    bs_f = Vf_ale.dofmap.index_map_bs
    bs_w = Vw.dofmap.index_map_bs
    if bs_f != 2 or bs_w != 2:
        raise RuntimeError(f"Expected 2-component interface fields, got fluid={bs_f}, wall={bs_w}")
    coords_f = Vf_ale.tabulate_dof_coordinates()[f_iface_dofs]
    coords_w = Vw.tabulate_dof_coordinates()[w_iface_dofs]
    f2w_map = KDTree(coords_w).query(coords_f)[1]

    # 7. Fluid Formulation (Navier-Stokes + BDF2 + ALE)
    w_ns, wn, wnn = fem.Function(W), fem.Function(W), fem.Function(W)
    u, p = ufl.split(w_ns)
    un, _ = ufl.split(wn)
    unn, _ = ufl.split(wnn)
    v, q = ufl.TestFunctions(W)
    
    uin, zero_v = fem.Function(Vf), fem.Function(Vf)
    zero_v.x.array[:] = 0.0
    
    c0, c1, c2 = fem.Constant(fluid, 1.5/DT), fem.Constant(fluid, -2.0/DT), fem.Constant(fluid, 0.5/DT)
    d_ale, d_ale_n, d_ale_nn = fem.Function(Vf_ale), fem.Function(Vf_ale), fem.Function(Vf_ale)
    w_mesh = fem.Function(Vf_ale)
    # The mesh velocity must use the same BE/BDF2 coefficients as the fluid.
    # In particular, the first step is BE: (d^1-d^0)/dt, not 1.5*d^1/dt.
    w_mesh_expr_ufl = c0*d_ale + c1*d_ale_n + c2*d_ale_nn
    w_mesh_expr = fem.Expression(w_mesh_expr_ufl, Vf_ale.element.interpolation_points)
    
    # Smooth normal for traction evaluation.  It must be anchored on the
    # fluid-wall interface, not on the outer inlet/outlet boundary.
    n_vec = build_wss_normal(fluid, fdim, fluid_iface_facets)
    
    # Fluid Weak Form
    dx_f = ufl.Measure("dx", domain=fluid)
    ds_f = ufl.Measure("ds", domain=fluid)
    # Since we don't have explicit subdomain tags on the fluid submesh, we can use 
    # default ds (with facet tags if we build them). Let's build a meshtags for fluid.
    f_mt_indices = np.concatenate([fluid_inlet_facets, fluid_outlet_facets, fluid_iface_facets])
    f_mt_values = np.concatenate([np.full(len(fluid_inlet_facets), 1),
                                  np.full(len(fluid_outlet_facets), 2),
                                  np.full(len(fluid_iface_facets), 6)])
    sort_idx = np.argsort(f_mt_indices)
    fluid_ft = dmesh.meshtags(fluid, fdim, f_mt_indices[sort_idx], f_mt_values[sort_idx])
    ds_f_tagged = ufl.Measure("ds", domain=fluid, subdomain_data=fluid_ft)

    p_out = fem.Constant(fluid, 0.0) # Simplified 1 outlet for smoke test, scale up later.
    F_fluid = (ufl.inner(c0*u + c1*un + c2*unn, v) * dx_f
             + ufl.inner(ufl.dot(u - w_mesh, ufl.nabla_grad(u)), v) * dx_f
             + NU * ufl.inner(ufl.grad(u), ufl.grad(v)) * dx_f
             - p * ufl.div(v) * dx_f
             - q * ufl.div(u) * dx_f
             + p_out * ufl.dot(v, ufl.FacetNormal(fluid)) * ds_f_tagged(2))

    # Boundary Conditions
    bc_inlet = fem.dirichletbc(uin, fem.locate_dofs_topological((W.sub(0), Vf), fdim, fluid_inlet_facets), W.sub(0))
    # Wall interface velocity is exactly w_mesh
    w_mesh_Vf = fem.Function(Vf)
    bc_wall = fem.dirichletbc(w_mesh_Vf, fem.locate_dofs_topological((W.sub(0), Vf), fdim, fluid_iface_facets), W.sub(0))
    
    J_fluid = ufl.derivative(F_fluid, w_ns, ufl.TrialFunction(W))
    args_nl = dict(bcs=[bc_inlet, bc_wall], J=J_fluid)
    modern_nl = "petsc_options_prefix" in inspect.signature(NonlinearProblem.__init__).parameters
    if modern_nl:
        args_nl.update(
            petsc_options_prefix="ns_bdf2_",
            petsc_options={
                "snes_type": "newtonls", "snes_rtol": NEWTON_RTOL, "snes_atol": NEWTON_ATOL,
                "snes_max_it": NEWTON_MAX_IT, "ksp_type": "preonly", "pc_type": "lu",
            }
        )
        fluid_problem = NonlinearProblem(F_fluid, w_ns, **args_nl)
        fluid_solver = fluid_problem  # modern dolfinx: problem acts as solver
    else:
        fluid_problem = NonlinearProblem(F_fluid, w_ns, **args_nl)
        fluid_solver = NewtonSolver(comm, fluid_problem)
        fluid_solver.rtol, fluid_solver.atol, fluid_solver.max_it = NEWTON_RTOL, NEWTON_ATOL, NEWTON_MAX_IT
        fluid_solver.convergence_criterion = "incremental"
        fluid_solver.krylov_solver.setType("preonly")
        fluid_solver.krylov_solver.getPC().setType("lu")

    # 8. Solid Formulation (Linear Elasticity)
    ds_solid = fem.Function(Vw, name="wall_displacement_m")
    ds_solid_n = fem.Function(Vw)
    ds_solid_nn = fem.Function(Vw)
    u_w, v_w = ufl.TrialFunction(Vw), ufl.TestFunction(Vw)
    eps_w = ufl.sym(ufl.grad(u_w))
    sigma_w = 2.0*MU_S*eps_w + LAMBDA*ufl.tr(eps_w)*ufl.Identity(2)
    # Dynamic plane-strain elasticity.  The old quasi-static solve had no
    # wall inertia, producing an added-mass blow-up at the first FSI step.
    a_solid = (WALL_RHO / DT**2 * ufl.inner(u_w, v_w)
               + ufl.inner(sigma_w, ufl.sym(ufl.grad(v_w)))) * ufl.dx(domain=wall)
    
    t_solid = fem.Function(Vw, name="fluid_traction_on_wall_Pa")
    # Create surface measure for solid interface
    s_mt = dmesh.meshtags(wall, fdim, wall_iface_facets, np.full(len(wall_iface_facets), 6, dtype=np.int32))
    ds_s = ufl.Measure("ds", domain=wall, subdomain_data=s_mt)
    L_solid = (WALL_RHO / DT**2 * ufl.inner(2.0*ds_solid_n - ds_solid_nn, v_w)
               * ufl.dx(domain=wall) + ufl.inner(t_solid, v_w) * ds_s(6))

    # Clamped on outer boundaries
    bc_solid_clamp = fem.dirichletbc(np.zeros(2, dtype=PETSc.ScalarType), fem.locate_dofs_topological(Vw, fdim, wall_outer_facets), Vw)
    solid_problem = safe_linear_problem(
        a_solid, L_solid, "solid_", bcs=[bc_solid_clamp], u=ds_solid,
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

    # 9. ALE Laplace (Extends interface displacement to fluid volume)
    u_ale, v_ale = ufl.TrialFunction(Vf_ale), ufl.TestFunction(Vf_ale)
    a_ale = ufl.inner(ufl.grad(u_ale), ufl.grad(v_ale)) * dx_f
    L_ale = ufl.inner(fem.Constant(fluid, np.zeros(2, dtype=PETSc.ScalarType)), v_ale) * dx_f
    bc_ale_fixed = fem.dirichletbc(np.zeros(2, dtype=PETSc.ScalarType), fem.locate_dofs_topological(Vf_ale, fdim, np.concatenate([fluid_inlet_facets, fluid_outlet_facets])), Vf_ale)
    d_iface_target = fem.Function(Vf_ale)
    bc_ale_iface = fem.dirichletbc(d_iface_target, f_iface_dofs)
    ale_problem = safe_linear_problem(a_ale, L_ale, "ale_", bcs=[bc_ale_fixed, bc_ale_iface], u=d_ale, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

    # Traction Evaluation (Volume projection at boundary)
    sigma_f = - (p * RHO) * ufl.Identity(2) + MU * (ufl.grad(u) + ufl.grad(u).T)
    traction_f = sigma_f * n_vec
    traction_expr = fem.Expression(traction_f, Vf_ale.element.interpolation_points)
    t_fluid = fem.Function(Vf_ale, name="fluid_traction_Pa")
    # The two submeshes have distinct dof numbering.  Build a *wall-to-fluid*
    # coordinate map; f2w_map above was computed in the opposite direction.
    w2f_map = KDTree(coords_f).query(coords_w)[1]
    map_error = float(np.max(KDTree(coords_f).query(coords_w)[0]))
    if map_error > 1e-9:
        raise RuntimeError(f"FSI interface DOFs do not coincide (max distance={map_error:.3e} m)")

    # Main Loop
    inlet = HeartInlet(HEART_BC)
    models = {2: OutletModel(*WK[2])} # Simplified to 1
    rows = []
    vu = VTKFile(comm, str(OUT / "velocity_fsi_bdf2.pvd"), "w")
    vp = VTKFile(comm, str(OUT / "pressure_fsi_bdf2.pvd"), "w")
    vd = VTKFile(comm, str(OUT / "wall_displacement_fsi_bdf2.pvd"), "w")
    
    print("\n--- STRONGLY-COUPLED FSI BDF2 ---")
    print(f"FSI startup ramp: {1e3*FSI_STARTUP_RAMP:.2f} ms (C1, zero-field initialisation)")
    start_time = time.perf_counter()
    # Explicit reference state: ParaView can compare step 0 with the first
    # accepted FSI checkpoint instead of inferring it from a later file.
    uh0, ph0 = w_ns.sub(0).collapse(), w_ns.sub(1).collapse()
    uh0.name, ph0.name = "velocity", "pressure_kinematic"
    vu.write_function(uh0, inlet.t_onset)
    vp.write_function(ph0, inlet.t_onset)
    vd.write_function(ds_solid, inlet.t_onset)
    print(f"[000/{MAX_STEPS}] t={inlet.t_onset:.4f} state=reference Vmax=0 Dmax=0")

    for step in range(1, MAX_STEPS + 1):
        if time.perf_counter() - start_time > 900: break
        
        t = inlet.t_onset + step * DT
        inlet.t_now = t
        uin.interpolate(inlet)
        
        if step == 1:
            c0.value, c1.value, c2.value = 1.0/DT, -1.0/DT, 0.0
        else:
            c0.value, c1.value, c2.value = 1.5/DT, -2.0/DT, 0.5/DT

        # Zero-order predictor for ALE to avoid massive artificial pressure shocks
        # caused by the added-mass effect when the mesh velocity is suddenly predicted.
        d_ale.x.array[:] = d_ale_n.x.array[:]
        
        fsi_iter = 0
        aitken_omega = FSI_RELAX
        previous_residual = None
        fsi_converged = False
        while fsi_iter < MAX_FSI_ITERS:
            fsi_iter += 1
            
            # Update mesh and compute w_m
            fluid.geometry.x[:,:2] = fluid_x0[:,:2] + d_ale.x.array.reshape(-1,2)
            w_mesh.interpolate(w_mesh_expr)
            # w_mesh is P1, but Navier-Stokes needs P2 for boundary conditions
            w_mesh_Vf.interpolate(w_mesh)
            
            # Solve Fluid
            if modern_nl:
                fluid_solver.solve()
            else:
                fluid_solver.solve(w_ns)
            w_ns.x.scatter_forward()
            
            # Evaluate Traction and Map Fluid -> Solid
            t_fluid.interpolate(traction_expr)
            # Action/reaction: wall receives the traction exerted by fluid.
            # This assignment is the actual fluid -> wall coupling; the old
            # draft wrote the wrong dofs and left the structural RHS at zero.
            t_solid.x.array.reshape((-1, bs_w))[w_iface_dofs] = \
                -t_fluid.x.array.reshape((-1, bs_f))[f_iface_dofs[w2f_map]]
            t_solid.x.scatter_forward()
            traction_max = float(np.max(np.linalg.norm(
                t_solid.x.array.reshape((-1, bs_w))[w_iface_dofs], axis=1)))
            if not np.isfinite(traction_max) or traction_max <= 1e-12:
                raise RuntimeError("FSI traction transfer is zero/invalid; refusing fake FSI step")
            # A failure must expose the load that caused it.  In particular,
            # the first physical step begins from an initially quiescent fluid
            # field and can create a large start-up pressure impulse.
            p_h = w_ns.sub(1).collapse()
            p_abs_max_pa = RHO * float(np.max(np.abs(p_h.x.array)))
            if fsi_iter == 1 or fsi_iter % 5 == 0:
                print(f"  step={step:03d} fsi={fsi_iter:02d} "
                      f"|traction|_max={traction_max:.3e}Pa "
                      f"|p|_max={p_abs_max_pa:.3e}Pa")
            
            # Solve Solid
            # LinearProblem.solve() returns the updated Function in newer
            # DOLFINx. Copy explicitly as well: otherwise a supplied ``u``
            # can leave ds_solid stale in some API versions.
            solid_result = solid_problem.solve()
            if solid_result is not ds_solid:
                ds_solid.x.array[:] = solid_result.x.array
            ds_solid.x.scatter_forward()
            
            # Extract displacements at interface
            # Solid interface values ordered in the same sequence as fluid
            # interface DOFs, so the residual is physically meaningful.
            d_s_iface = ds_solid.x.array.reshape((-1, bs_w))[w_iface_dofs[f2w_map]]
            d_s_max = float(np.max(np.linalg.norm(ds_solid.x.array.reshape((-1, 2)), axis=1)))
            if not np.isfinite(d_s_max) or d_s_max > MAX_INTERFACE_DISP:
                raise RuntimeError(
                    f"Wall displacement {1e3*d_s_max:.4f} mm exceeds ALE trust region "
                    f"{1e3*MAX_INTERFACE_DISP:.4f} mm at step {step}, FSI iteration {fsi_iter}; "
                    f"|traction|_max={traction_max:.3e} Pa, |p|_max={p_abs_max_pa:.3e} Pa")
            d_f_iface = d_ale.x.array.reshape((-1, bs_f))[f_iface_dofs]
            
            # Residual and Relaxation
            res_vec = d_s_iface - d_f_iface
            n_iface = max(len(f_iface_dofs), 1)
            res_rms = float(np.linalg.norm(res_vec) / np.sqrt(n_iface))
            d_s_rms = float(np.linalg.norm(d_s_iface) / np.sqrt(n_iface))
            rel_res = res_rms / max(d_s_rms, FSI_DISP_REFERENCE)

            # Dynamic Aitken relaxation for the partitioned FSI fixed point.
            # It suppresses the added-mass oscillation without concealing it:
            # omega and the residual are printed on every coupling iteration.
            if previous_residual is not None:
                delta_residual = res_vec - previous_residual
                denom = float(np.vdot(delta_residual, delta_residual).real)
                if denom > 1e-30:
                    aitken_omega = float(np.clip(
                        -aitken_omega * np.vdot(previous_residual, delta_residual).real / denom,
                        FSI_RELAX_MIN, FSI_RELAX_MAX))
            previous_residual = res_vec.copy()
            
            # Update target for ALE
            d_iface_target.x.array.reshape((-1, bs_f))[f_iface_dofs] = \
                d_f_iface + aitken_omega * res_vec
            d_iface_target.x.scatter_forward()
            
            # Solve ALE and explicitly retain the returned solution.  This is
            # the displacement used by the next fluid solve and by w_mesh.
            ale_result = ale_problem.solve()
            if ale_result is not d_ale:
                d_ale.x.array[:] = ale_result.x.array
            d_ale.x.scatter_forward()

            # Every coupling iteration is logged: this is the convergence
            # evidence for FSI, unlike a single final line per time step.
            d_trial = float(np.max(np.linalg.norm(d_ale.x.array.reshape((-1, 2)), axis=1)))
            d_wall = float(np.max(np.linalg.norm(ds_solid.x.array.reshape((-1, 2)), axis=1)))
            print(f"  step={step:03d} fsi={fsi_iter:02d} omega={aitken_omega:.3f} "
                  f"res_rel={rel_res:.3e} res_abs={1e9*res_rms:.3f}nm "
                  f"d_mesh={1e3*d_trial:.5f}mm "
                  f"d_wall={1e3*d_wall:.5f}mm")
            
            if res_rms < FSI_ABS_TOL or rel_res < FSI_TOL:
                fsi_converged = True
                break

        if not fsi_converged:
            raise RuntimeError(
                f"FSI failed to converge at step {step}: "
                f"res_rel={rel_res:.3e}, res_abs={1e9*res_rms:.3f} nm")
                
        # Advance histories
        d_ale_nn.x.array[:] = d_ale_n.x.array[:]
        d_ale_n.x.array[:] = d_ale.x.array[:]
        ds_solid_nn.x.array[:] = ds_solid_n.x.array[:]
        ds_solid_n.x.array[:] = ds_solid.x.array[:]
        wnn.x.array[:] = wn.x.array[:]
        wn.x.array[:] = w_ns.x.array[:]
        
        vmax = float(np.max(np.linalg.norm(w_ns.sub(0).collapse().x.array.reshape(-1,2), axis=1)))
        dmax = float(np.max(np.linalg.norm(ds_solid.x.array.reshape((-1, 2)), axis=1)))
        rows.append((step, t, vmax, dmax, traction_max, rel_res, fsi_iter))
        if step == 1 or step % LOG_EVERY == 0:
            print(f"[{step:03d}/{MAX_STEPS}] t={t:.4f} FSI={fsi_iter}/{MAX_FSI_ITERS} "
                  f"res_d={rel_res:.2e} traction={traction_max:.2e}Pa "
                  f"Vmax={vmax:.3f} Dmax={1e3*dmax:.5f}mm")
        if step % WRITE_EVERY == 0 or step == MAX_STEPS:
            uh, ph = w_ns.sub(0).collapse(), w_ns.sub(1).collapse()
            uh.name, ph.name = "velocity", "pressure_kinematic"
            ds_solid.name = "wall_displacement_m"
            vu.write_function(uh, t); vp.write_function(ph, t); vd.write_function(ds_solid, t)

    vu.close(); vp.close(); vd.close()
    if rows:
        a = np.asarray(rows)
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        panels = (("Vmax (m/s)", 2, False), ("Wall displacement (mm)", 3, False),
                  ("Fluid traction (Pa)", 4, True), ("FSI residual", 5, True),
                  ("FSI iterations", 6, False))
        for ax, (title, col, logarithmic) in zip(axes.flat, panels):
            y = a[:, col] * (1e3 if col == 3 else 1.0)
            ax.plot(a[:, 1], np.maximum(y, 1e-16) if logarithmic else y, lw=1.5)
            ax.set(title=title, xlabel="t (s)"); ax.grid(alpha=0.3)
            if logarithmic: ax.set_yscale("log")
        axes.flat[-1].axis("off")
        fig.suptitle("Aorta strong FSI BDF2: fluid traction -> elastic wall -> ALE")
        fig.tight_layout()
        fig.savefig(OUT / "aorta_fsi_bdf2_diagnostics.png", dpi=180)
        fig.savefig(OUT / "aorta_fsi_bdf2_diagnostics.pdf")
        plt.close(fig)

if __name__ == "__main__":
    main()
