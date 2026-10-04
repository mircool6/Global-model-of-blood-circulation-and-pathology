import gmsh
import numpy as np
import os
import matplotlib.pyplot as plt

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem
from dolfinx.io import gmsh as gmshio, VTKFile
from basix.ufl import element, mixed_element
import ufl
from ufl import TrialFunction, TestFunctions, split, inner, grad, div, dot, cross

# ==========================================================
# PHYSICAL CONSTANTS (mm-kg-s system)
# ==========================================================
rho = 1.05e-6  # kg/mm^3
mu  = 3.5e-6   # kg/(mm*s)
nu  = PETSc.ScalarType(mu / rho)

STEP_PATH   = "/mnt/c/Users/Mirco/Downloads/turbine_good_scketch1.step"
TEMP_MESH   = "/mnt/c/Users/Mirco/Downloads/temp_sweep_mesh.msh"
OUTPUT_DIR  = "/mnt/c/Users/Mirco/Downloads"

# Stagger mesh path (pre-built mesh)
STAGGER_MESH = "/mnt/c/Work/Cursor_CFD/salome_mesh/stagger1mesh.msh"

# Original baseline values from calculate_metrics.py run
BASELINE_HEMOLYSIS  = 207.1  # W
BASELINE_STAGNATION = 0.19   # %
BASELINE_TORQUE     = 0.48   # N*mm


def make_wedge_tool(a_start, a_end, R=10.0, H=55.0):
    """
    Creates a wedge-shaped solid covering the angular sector
    from a_start to a_end (in radians) along the Z axis.
    """
    pts = [
        gmsh.model.occ.addPoint(0, 0, -1),
        gmsh.model.occ.addPoint(R * np.cos(a_start), R * np.sin(a_start), -1),
        gmsh.model.occ.addPoint(R * np.cos(a_end),   R * np.sin(a_end),   -1),
        gmsh.model.occ.addPoint(0, 0, H),
        gmsh.model.occ.addPoint(R * np.cos(a_start), R * np.sin(a_start), H),
        gmsh.model.occ.addPoint(R * np.cos(a_end),   R * np.sin(a_end),   H),
    ]
    ls = [
        gmsh.model.occ.addLine(pts[0], pts[1]),  # 0
        gmsh.model.occ.addLine(pts[1], pts[2]),  # 1
        gmsh.model.occ.addLine(pts[2], pts[0]),  # 2
        gmsh.model.occ.addLine(pts[3], pts[4]),  # 3
        gmsh.model.occ.addLine(pts[4], pts[5]),  # 4
        gmsh.model.occ.addLine(pts[5], pts[3]),  # 5
        gmsh.model.occ.addLine(pts[0], pts[3]),  # 6
        gmsh.model.occ.addLine(pts[1], pts[4]),  # 7
        gmsh.model.occ.addLine(pts[2], pts[5]),  # 8
    ]
    cls = [
        gmsh.model.occ.addCurveLoop([ls[0], ls[1], ls[2]]),          # bottom
        gmsh.model.occ.addCurveLoop([ls[3], ls[4], ls[5]]),          # top
        gmsh.model.occ.addCurveLoop([ls[0], ls[7], -ls[3], -ls[6]]),  # face 1
        gmsh.model.occ.addCurveLoop([ls[1], ls[8], -ls[4], -ls[7]]),  # face 2
        gmsh.model.occ.addCurveLoop([ls[2], ls[6], -ls[5], -ls[8]]),  # face 3
    ]
    surfs = [gmsh.model.occ.addPlaneSurface([cl]) for cl in cls]
    shell = gmsh.model.occ.addSurfaceLoop(surfs)
    return gmsh.model.occ.addVolume([shell])


def build_fluid_mesh_from_step(n_blades, mesh_path, mesh_min=0.3, mesh_max=1.6):
    """
    Imports the original STEP, keeps only n_blades out of 3
    by intersecting with angular wedge sectors, then subtracts from
    the fluid cylinder and generates the tetrahedral mesh.
    """
    print(f"\n[GMSH] Building fluid mesh from original STEP with N_BLADES = {n_blades}...")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("fluid_sweep")

    # 1. Import original STEP
    shapes = gmsh.model.occ.importShapes(STEP_PATH)
    gmsh.model.occ.synchronize()

    # 2. Identify casing vs rotor
    rotor_tag = None
    casing_tag = None
    for dim, tag in gmsh.model.getEntities(3):
        bbox = gmsh.model.getBoundingBox(dim, tag)
        if bbox[5] < 45.0:
            rotor_tag = tag
        else:
            casing_tag = tag

    print(f"  Rotor tag: {rotor_tag}, Casing tag: {casing_tag}")

    # 3. Scale rotor 0.8 in XY (same as generate_meshes.py)
    gmsh.model.occ.dilate([(3, rotor_tag)], 0, 0, 0, 0.8, 0.8, 1.0)
    gmsh.model.occ.synchronize()

    # 4. If fewer than 3 blades, intersect with union of N wedge sectors
    HUB_RADIUS = 1.5  # mm — shaft/hub radius (blades start at ~1.5 mm)

    if n_blades < 3:
        print(f"  Trimming to {n_blades} blade(s), keeping hub intact...")

        # Step A: Extract hub separately (cylinder R=HUB_RADIUS, full length)
        hub_cyl = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 42.0, HUB_RADIUS)
        gmsh.model.occ.synchronize()

        # Intersect rotor with hub cylinder → get pure hub solid
        hub_result, _ = gmsh.model.occ.intersect(
            [(3, rotor_tag)], [(3, hub_cyl)],
            removeObject=False, removeTool=True
        )
        gmsh.model.occ.synchronize()
        hub_tag = hub_result[0][1]
        hub_mass = gmsh.model.occ.getMass(3, hub_tag)
        print(f"  Hub extracted: volume = {hub_mass:.2f} mm^3")

        # Step B: Cut hub out of rotor to get blades-only solid
        blades_result, _ = gmsh.model.occ.cut(
            [(3, rotor_tag)], [(3, hub_tag)],
            removeObject=True, removeTool=False
        )
        gmsh.model.occ.synchronize()
        blades_tag = blades_result[0][1]
        blades_mass = gmsh.model.occ.getMass(3, blades_tag)
        print(f"  Blades-only solid: volume = {blades_mass:.2f} mm^3")

        # Step C: Intersect blades-only with N wedge sectors
        blade_sectors = [
            (np.radians(-67.0),  np.radians(53.0)),
            (np.radians(173.0),  np.radians(293.0)),
            (np.radians(53.0),   np.radians(173.0)),
        ]
        sector_tags = []
        for i in range(n_blades):
            a_start, a_end = blade_sectors[i]
            v = make_wedge_tool(a_start, a_end)
            sector_tags.append((3, v))

        gmsh.model.occ.synchronize()

        if len(sector_tags) > 1:
            united, _ = gmsh.model.occ.fuse([sector_tags[0]], sector_tags[1:])
            tool = united
        else:
            tool = sector_tags

        gmsh.model.occ.synchronize()

        selected_blades, _ = gmsh.model.occ.intersect(
            [(3, blades_tag)], tool,
            removeObject=True, removeTool=True
        )
        gmsh.model.occ.synchronize()
        sel_mass = gmsh.model.occ.getMass(selected_blades[0][0], selected_blades[0][1])
        print(f"  Selected {n_blades} blade(s): volume = {sel_mass:.2f} mm^3")

        # Step D: Fuse selected blades back with hub
        fused_rotor, _ = gmsh.model.occ.fuse(
            [(3, hub_tag)], selected_blades,
            removeObject=True, removeTool=True
        )
        gmsh.model.occ.synchronize()
        rotor_tag = fused_rotor[0][1]
        final_mass = gmsh.model.occ.getMass(3, rotor_tag)
        print(f"  Final rotor (hub + {n_blades} blade(s)): volume = {final_mass:.2f} mm^3")

    else:
        mass = gmsh.model.occ.getMass(3, rotor_tag)
        print(f"  Full rotor kept (3 blades), volume: {mass:.2f} mm^3")

    # 5. Create the fluid cylinder casing (R=4.5 mm, H=50 mm)
    fluid_cyl = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 50.0, 4.5)
    gmsh.model.occ.synchronize()

    # 6. Boolean subtraction: fluid = cylinder - rotor
    fluid_entities, _ = gmsh.model.occ.cut(
        [(3, fluid_cyl)], [(3, rotor_tag)],
        removeTool=True, removeObject=True
    )
    gmsh.model.occ.synchronize()

    if not fluid_entities:
        raise RuntimeError("Boolean cut produced no fluid volume!")

    fluid_tag = fluid_entities[0][1]
    print(f"  Fluid domain tag: {fluid_tag}")

    # 7. Classify boundary surfaces geometrically
    boundary = gmsh.model.getBoundary(
        [(3, fluid_tag)], combined=False, oriented=False, recursive=False
    )

    inlet_surfaces    = []
    outlet_surfaces   = []
    wall_surfaces     = []
    turbine_surfaces  = []

    for dim, tag in boundary:
        bbox = gmsh.model.getBoundingBox(dim, tag)
        zmin, zmax = bbox[2], bbox[5]
        xmin, xmax = bbox[0], bbox[3]
        ymin, ymax = bbox[1], bbox[4]
        max_r = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))

        if abs(zmax - 50.0) < 0.5 and (zmax - zmin) < 0.5:
            inlet_surfaces.append(tag)
        elif abs(zmin - 0.0) < 0.5 and (zmax - zmin) < 0.5:
            if max_r > 3.5:
                outlet_surfaces.append(tag)
            else:
                turbine_surfaces.append(tag)
        else:
            if max_r > 4.0:
                wall_surfaces.append(tag)
            else:
                turbine_surfaces.append(tag)

    print(f"  Inlet: {len(inlet_surfaces)}, Outlet: {len(outlet_surfaces)}, "
          f"Wall: {len(wall_surfaces)}, Turbine: {len(turbine_surfaces)}")

    # 8. Physical groups
    gmsh.model.addPhysicalGroup(2, inlet_surfaces,   tag=325, name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfaces,  tag=324, name="outlet")
    gmsh.model.addPhysicalGroup(2, wall_surfaces,    tag=326, name="wall")
    gmsh.model.addPhysicalGroup(2, turbine_surfaces, tag=327, name="turbine")
    gmsh.model.addPhysicalGroup(3, [fluid_tag],      tag=328, name="fluid")

    # 9. Sizing field: refined near turbine surfaces
    if turbine_surfaces:
        dist_field = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(dist_field, "SurfacesList", turbine_surfaces)
        gmsh.model.mesh.field.setNumber(dist_field, "Sampling", 100)

        thresh_field = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(thresh_field, "InField",  dist_field)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMin",  mesh_min)
        gmsh.model.mesh.field.setNumber(thresh_field, "SizeMax",  mesh_max)
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMin",  0.5)
        gmsh.model.mesh.field.setNumber(thresh_field, "DistMax",  3.0)
        gmsh.model.mesh.field.setAsBackgroundMesh(thresh_field)

    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_min)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_max)
    factor = mesh_max / 1.6
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", factor)

    # 10. Mesh
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write(mesh_path)
    gmsh.finalize()
    print(f"[GMSH] Mesh saved to {mesh_path}")


def solve_cfd_and_get_metrics(mesh_path, n_blades):
    """
    Loads mesh, solves Stokes + NS, saves VTU fields, returns metrics.
    """
    print("[FENICS] Solving Stokes + Navier-Stokes...")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(mesh_path)
    mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
    gmsh.finalize()

    mesh      = mesh_data.mesh
    facet_tags = mesh_data.facet_tags
    mesh.topology.create_connectivity(mesh.topology.dim - 1, 0)

    cell = mesh.topology.cell_name()
    V_el = element("Lagrange", cell, 2, shape=(3,))
    Q_el = element("Lagrange", cell, 1)
    W    = fem.functionspace(mesh, mixed_element([V_el, Q_el]))
    print(f"  Mixed DOFs: {W.dofmap.index_map.size_global}")

    dx_m = ufl.Measure("dx", domain=mesh)
    ds_m = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

    # ---- Boundary Conditions ----
    def vector_bc(W_sub, facets, val):
        V_c, _ = W_sub.collapse()
        u_bc   = fem.Function(V_c)
        u_bc.interpolate(lambda x: np.tile(val, (x.shape[1], 1)).T)
        dofs = fem.locate_dofs_topological((W_sub, V_c), mesh.topology.dim - 1, facets)
        return fem.dirichletbc(u_bc, dofs, W_sub)

    Uin  = 100.0  # mm/s
    bcs  = []
    bcs.append(vector_bc(W.sub(0), facet_tags.find(325),
                         np.array([0.0, 0.0, -Uin], dtype=PETSc.ScalarType)))
    wall_f = np.concatenate([facet_tags.find(326), facet_tags.find(327)])
    bcs.append(vector_bc(W.sub(0), wall_f,
                         np.zeros(3, dtype=PETSc.ScalarType)))
    p_dofs = fem.locate_dofs_topological(W.sub(1), mesh.topology.dim - 1,
                                         facet_tags.find(324))
    bcs.append(fem.dirichletbc(PETSc.ScalarType(0), p_dofs, W.sub(1)))

    # ---- Stokes solve ----
    w_s  = TrialFunction(W)
    v, q = TestFunctions(W)
    u_s, p_s = split(w_s)

    a = (nu * inner(grad(u_s), grad(v)) - p_s * div(v) - q * div(u_s)) * dx_m
    L = (inner(fem.Constant(mesh, PETSc.ScalarType((0, 0, 0))), v)
         + fem.Constant(mesh, PETSc.ScalarType(0.0)) * q) * dx_m

    from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
    stokes = LinearProblem(a, L, bcs=bcs,
                           petsc_options_prefix="stokes_sol_",
                           petsc_options={
                               "ksp_type": "preonly",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps",
                               "stokes_sol_mat_mumps_icntl_14": 120,
                               "stokes_sol_mat_mumps_icntl_23": 2000
                           })
    w_st = stokes.solve()
    print("  Stokes solved.")

    # ---- Navier-Stokes solve ----
    w_ns   = fem.Function(W)
    w_ns.x.array[:] = w_st.x.array[:]
    u_ns, p_ns = split(w_ns)

    F = (inner(dot(grad(u_ns), u_ns), v)
         + nu * inner(grad(u_ns), grad(v))
         - p_ns * div(v)
         - q * div(u_ns)) * dx_m
    J = ufl.derivative(F, w_ns)

    ns = NonlinearProblem(F, w_ns, J=J, bcs=bcs,
                          petsc_options_prefix="ns_sol_",
                          petsc_options={
                              "snes_type": "newtonls",
                              "snes_linesearch_type": "basic",
                              "snes_max_it": 50,
                              "snes_rtol": 1e-4,
                              "snes_atol": 1e-4,
                              "snes_monitor": None,
                              "snes_converged_reason": None,
                              "ksp_type": "preonly",
                              "pc_type": "lu",
                              "pc_factor_mat_solver_type": "mumps",
                              "ns_sol_mat_mumps_icntl_14": 120,
                              "ns_sol_mat_mumps_icntl_23": 2000
                          })
    w_ns_sol = ns.solve()
    print("  Navier-Stokes solved.")

    # ---- Split and compute pressure ----
    u_sol, phi_sol = w_ns_sol.split()
    u_sol   = u_sol.collapse()
    phi_sol = phi_sol.collapse()
    p_sol   = fem.Function(phi_sol.function_space, name="pressure")
    p_sol.x.array[:] = phi_sol.x.array[:] * rho

    # ---- Save VTU for ParaView ----
    u_sol.name = "velocity"
    vel_path = f"{OUTPUT_DIR}/velocity_sweep_N{n_blades}.pvd"
    pre_path = f"{OUTPUT_DIR}/pressure_sweep_N{n_blades}.pvd"
    with VTKFile(mesh.comm, vel_path, "w") as f:
        f.write_function(u_sol, 0.0)
    with VTKFile(mesh.comm, pre_path, "w") as f:
        f.write_function(p_sol, 0.0)
    print(f"  VTU saved: {vel_path}")

    # ---- Metrics ----
    vol_total = fem.assemble_scalar(
        fem.form(fem.Constant(mesh, PETSc.ScalarType(1.0)) * dx_m))

    S      = 0.5 * (grad(u_sol) + grad(u_sol).T)
    e_diss = fem.assemble_scalar(
        fem.form(2.0 * fem.Constant(mesh, PETSc.ScalarType(mu)) * inner(S, S) * dx_m))

    u_mag  = ufl.sqrt(inner(u_sol, u_sol))
    v_stag = fem.assemble_scalar(
        fem.form(ufl.conditional(u_mag < 10.0, 1.0, 0.0) * dx_m))
    stag_frac = (v_stag / vol_total) * 100.0

    n_fac  = ufl.FacetNormal(mesh)
    sigma  = -p_sol * ufl.Identity(3) + 2.0 * fem.Constant(mesh, PETSc.ScalarType(mu)) * S
    x_c    = ufl.SpatialCoordinate(mesh)
    r_vec  = ufl.as_vector([x_c[0], x_c[1], 0.0])
    T_z    = fem.assemble_scalar(
        fem.form(dot(cross(r_vec, dot(sigma, n_fac)),
                     ufl.as_vector([0.0, 0.0, 1.0])) * ds_m(327)))

    n_dofs = W.dofmap.index_map.size_global

    print(f"  [RESULT] DOFs:       {n_dofs}")
    print(f"  [RESULT] Hemolysis:  {e_diss:.3e} W")
    print(f"  [RESULT] Stagnation: {stag_frac:.2f} %")
    print(f"  [RESULT] Torque:     {T_z:.3e} N*mm")

    return e_diss, stag_frac, T_z, n_dofs


def solve_cfd_on_stagger_mesh():
    """
    Loads pre-built stagger mesh (stagger1mesh.msh), solves Stokes + NS,
    saves VTU fields, returns metrics.
    """
    print("[FENICS] Solving Stokes + Navier-Stokes on STAGGER mesh...")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(STAGGER_MESH)
    
    # Dynamically extract physical tags by name
    groups = gmsh.model.getPhysicalGroups()
    tag_map = {}
    for dim, tag in groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        if name:
            tag_map[name] = tag
            
    print(f"  Physical groups found: {tag_map}")
    
    mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
    gmsh.finalize()

    mesh      = mesh_data.mesh
    facet_tags = mesh_data.facet_tags
    mesh.topology.create_connectivity(mesh.topology.dim - 1, 0)

    cell = mesh.topology.cell_name()
    V_el = element("Lagrange", cell, 2, shape=(3,))
    Q_el = element("Lagrange", cell, 1)
    W    = fem.functionspace(mesh, mixed_element([V_el, Q_el]))
    print(f"  Mixed DOFs: {W.dofmap.index_map.size_global}")

    # Resolve actual tags
    tag_inlet = tag_map.get("inlet", 326)
    tag_outlet = tag_map.get("outlet", 324)
    tag_wall = tag_map.get("wall", 326)
    tag_turbine = tag_map.get("turbine", 327)

    dx_m = ufl.Measure("dx", domain=mesh)
    ds_m = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

    # ---- Boundary Conditions ----
    def vector_bc(W_sub, facets, val):
        V_c, _ = W_sub.collapse()
        u_bc   = fem.Function(V_c)
        u_bc.interpolate(lambda x: np.tile(val, (x.shape[1], 1)).T)
        dofs = fem.locate_dofs_topological((W_sub, V_c), mesh.topology.dim - 1, facets)
        return fem.dirichletbc(u_bc, dofs, W_sub)

    Uin  = 100.0  # mm/s
    bcs  = []
    bcs.append(vector_bc(W.sub(0), facet_tags.find(tag_inlet),
                         np.array([0.0, 0.0, -Uin], dtype=PETSc.ScalarType)))
    
    wall_f = np.concatenate([facet_tags.find(tag_wall), facet_tags.find(tag_turbine)])
    bcs.append(vector_bc(W.sub(0), wall_f,
                         np.zeros(3, dtype=PETSc.ScalarType)))
                        
    p_dofs = fem.locate_dofs_topological(W.sub(1), mesh.topology.dim - 1,
                                         facet_tags.find(tag_outlet))
    bcs.append(fem.dirichletbc(PETSc.ScalarType(0), p_dofs, W.sub(1)))

    # ---- Stokes solve ----
    w_s  = TrialFunction(W)
    v, q = TestFunctions(W)
    u_s, p_s = split(w_s)

    a = (nu * inner(grad(u_s), grad(v)) - p_s * div(v) - q * div(u_s)) * dx_m
    L = (inner(fem.Constant(mesh, PETSc.ScalarType((0, 0, 0))), v)
         + fem.Constant(mesh, PETSc.ScalarType(0.0)) * q) * dx_m

    from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
    stokes_problem = LinearProblem(a, L, bcs=bcs,
                           petsc_options_prefix="stokes_sol_",
                           petsc_options={
                               "ksp_type": "preonly",
                               "pc_type": "lu",
                               "pc_factor_mat_solver_type": "mumps",
                               "stokes_sol_mat_mumps_icntl_14": 120,
                               "stokes_sol_mat_mumps_icntl_23": 2000
                           })
    w_st = stokes_problem.solve()
    print("  Stokes solved.")

    # ---- Navier-Stokes solve ----
    w_ns   = fem.Function(W)
    w_ns.x.array[:] = w_st.x.array[:]
    u_ns, p_ns = split(w_ns)

    F = (inner(dot(grad(u_ns), u_ns), v)
         + nu * inner(grad(u_ns), grad(v))
         - p_ns * div(v)
         - q * div(u_ns)) * dx_m
    J = ufl.derivative(F, w_ns)

    ns = NonlinearProblem(F, w_ns, J=J, bcs=bcs,
                          petsc_options_prefix="ns_sol_",
                          petsc_options={
                              "snes_type": "newtonls",
                              "snes_linesearch_type": "basic",
                              "snes_max_it": 50,
                              "snes_rtol": 1e-4,
                              "snes_atol": 1e-4,
                              "ksp_type": "preonly",
                              "pc_type": "lu",
                              "pc_factor_mat_solver_type": "mumps",
                              "ns_sol_mat_mumps_icntl_14": 120,
                              "ns_sol_mat_mumps_icntl_23": 2000
                          })
    w_ns_sol = ns.solve()
    print("  Navier-Stokes solved.")

    # ---- Split and compute pressure ----
    u_sol, phi_sol = w_ns_sol.split()
    u_sol   = u_sol.collapse()
    phi_sol = phi_sol.collapse()
    p_sol   = fem.Function(phi_sol.function_space, name="pressure")
    p_sol.x.array[:] = phi_sol.x.array[:] * rho

    # ---- Save VTU for ParaView ----
    u_sol.name = "velocity"
    vel_path = f"{OUTPUT_DIR}/velocity_stagger.pvd"
    pre_path = f"{OUTPUT_DIR}/pressure_stagger.pvd"
    with VTKFile(mesh.comm, vel_path, "w") as f:
        f.write_function(u_sol, 0.0)
    with VTKFile(mesh.comm, pre_path, "w") as f:
        f.write_function(p_sol, 0.0)
    print(f"  VTU saved: {vel_path}")

    # ---- Metrics ----
    vol_total = fem.assemble_scalar(
        fem.form(fem.Constant(mesh, PETSc.ScalarType(1.0)) * dx_m))

    S      = 0.5 * (grad(u_sol) + grad(u_sol).T)
    e_diss = fem.assemble_scalar(
        fem.form(2.0 * fem.Constant(mesh, PETSc.ScalarType(mu)) * inner(S, S) * dx_m))

    u_mag  = ufl.sqrt(inner(u_sol, u_sol))
    v_stag = fem.assemble_scalar(
        fem.form(ufl.conditional(u_mag < 10.0, 1.0, 0.0) * dx_m))
    stag_frac = (v_stag / vol_total) * 100.0

    n_fac  = ufl.FacetNormal(mesh)
    sigma  = -p_sol * ufl.Identity(3) + 2.0 * fem.Constant(mesh, PETSc.ScalarType(mu)) * S
    x_c    = ufl.SpatialCoordinate(mesh)
    r_vec  = ufl.as_vector([x_c[0], x_c[1], 0.0])
    T_z    = fem.assemble_scalar(
        fem.form(dot(cross(r_vec, dot(sigma, n_fac)),
                     ufl.as_vector([0.0, 0.0, 1.0])) * ds_m(tag_turbine)))

    n_dofs = W.dofmap.index_map.size_global

    print(f"  [RESULT] DOFs:       {n_dofs}")
    print(f"  [RESULT] Hemolysis:  {e_diss:.3e} W")
    print(f"  [RESULT] Stagnation: {stag_frac:.2f} %")
    print(f"  [RESULT] Torque:     {T_z:.3e} N*mm")

    return e_diss, stag_frac, T_z, n_dofs


def run_sweep():
    # N=1 (min), N=2 (mid), N=3 (original full rotor)
    blade_counts    = [1, 2, 3]
    hemolysis_list  = []
    stagnation_list = []
    torque_list     = []
    dof_list        = []

    # To fairly compare blade geometries, all N run with IDENTICAL mesh sizing.
    _h_min, _h_max = 0.42, 1.90   # uniform for all blade counts
    sweep_mesh_params = {
        1: (_h_min, _h_max),
        2: (_h_min, _h_max),
        3: (_h_min, _h_max),
    }

    for n in blade_counts:
        mesh_min, mesh_max = sweep_mesh_params[n]
        build_fluid_mesh_from_step(n, TEMP_MESH,
                                   mesh_min=mesh_min,
                                   mesh_max=mesh_max)
        h, s, t, d = solve_cfd_and_get_metrics(TEMP_MESH, n)
        hemolysis_list.append(h)
        stagnation_list.append(s)
        torque_list.append(abs(t))
        dof_list.append(d)

    if os.path.exists(TEMP_MESH):
        os.remove(TEMP_MESH)

    # Add stagger mesh results (pre-built mesh) - separate point
    print("\n[SOLVING] Stagger mesh (pre-built from NS_Gimini_stagger)...")
    h_stag, s_stag, t_stag, d_stag = solve_cfd_on_stagger_mesh()
    
    # Add stagger as separate point (after 3 blades, with its DOF)
    blade_counts.append(3.5)  # Use 3.5 to place it after 3 blades
    hemolysis_list.append(h_stag)
    stagnation_list.append(s_stag)
    torque_list.append(abs(t_stag))
    dof_list.append(d_stag)

    print("\n==========================================================")
    print("SWEEP COMPLETE — GENERATING PLOTS")
    print("==========================================================")

    # ---- Plot 1: Safety Metrics ----
    fig, ax1 = plt.subplots(figsize=(9, 6))
    c1 = 'tab:red'
    ax1.set_xlabel('Number of Blades', fontweight='bold', fontsize=12)
    ax1.set_ylabel('Hemolysis: Dissipation Rate (W)', color=c1,
                   fontweight='bold', fontsize=12)
    l1 = ax1.plot(blade_counts, hemolysis_list, 'o-',
                  color=c1, linewidth=2.5, markersize=9, label='Hemolysis (parametric)')
    ax1.axhline(BASELINE_HEMOLYSIS, color=c1, linestyle='--', linewidth=1.5,
                label=f'Original 3-blade ({BASELINE_HEMOLYSIS} W)')
    ax1.tick_params(axis='y', labelcolor=c1)
    
    # Stagger mesh as separate point
    ax1.plot(3.5, h_stag, 's', color='tab:orange', markersize=12,
             label=f'Stagger mesh ({d_stag:,} DOF)')
    
    ax1.set_xticks(blade_counts)
    ax1.set_xticklabels(['1 blade\n(min)', '2 blades\n(mid)', '3 blades\n(original)', 'Stagger\nmesh'])

    ax2 = ax1.twinx()
    c2 = 'tab:blue'
    ax2.set_ylabel('Thrombosis: Stagnant Volume (%)', color=c2,
                   fontweight='bold', fontsize=12)
    l2 = ax2.plot(blade_counts, stagnation_list, 's-',
                  color=c2, linewidth=2.5, markersize=9, label='Stagnation (parametric)')
    ax2.axhline(BASELINE_STAGNATION, color=c2, linestyle='--', linewidth=1.5,
                label=f'Original 3-blade ({BASELINE_STAGNATION}%)')
    ax2.tick_params(axis='y', labelcolor=c2)
    
    # Stagger mesh as separate point
    ax2.plot(3.5, s_stag, 's', color='tab:cyan', markersize=12)

    handles = l1 + [plt.Line2D([0], [0], color=c1, linestyle='--', label=f'Original baseline ({BASELINE_HEMOLYSIS} W)')] + \
              [plt.Line2D([0], [0], marker='s', color='tab:orange', linestyle='None', markersize=10, label=f'Stagger mesh ({d_stag:,} DOF)')] + \
              l2 + [plt.Line2D([0], [0], color=c2, linestyle='--', label=f'Original baseline ({BASELINE_STAGNATION}%)')]
    ax1.legend(handles, [h.get_label() for h in handles], loc='upper left', fontsize=9)

    plt.title('Blood Pump Safety: Blade Count vs Hemolysis & Thrombosis',
              fontweight='bold', fontsize=13, pad=15)
    plt.tight_layout()
    p1 = f"{OUTPUT_DIR}/blade_sweep_safety.png"
    plt.savefig(p1, dpi=300)
    print(f"Saved: {p1}")

    # ---- Plot 2: Torque ----
    fig2, ax = plt.subplots(figsize=(9, 6))
    ax.plot(blade_counts, torque_list, '^-',
            color='tab:green', linewidth=2.5, markersize=9, label='Torque (parametric)')
    ax.axhline(BASELINE_TORQUE, color='tab:green', linestyle='--', linewidth=1.5,
               label=f'Original 3-blade ({BASELINE_TORQUE} N*mm)')
    ax.plot(3.5, abs(t_stag), 's', color='tab:purple', markersize=12,
            label=f'Stagger mesh ({d_stag:,} DOF)')
    ax.set_xlabel('Number of Blades', fontweight='bold', fontsize=12)
    ax.set_ylabel('Torque Magnitude (N·mm)', fontweight='bold', fontsize=12)
    ax.set_xticks(blade_counts)
    ax.set_xticklabels(['1 blade\n(min)', '2 blades\n(mid)', '3 blades\n(original)', 'Stagger\nmesh'])
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.title('Blood Pump Efficiency: Blade Count vs Torque',
              fontweight='bold', fontsize=13, pad=15)
    plt.tight_layout()
    p2 = f"{OUTPUT_DIR}/blade_sweep_torque.png"
    plt.savefig(p2, dpi=300)
    print(f"Saved: {p2}")


def mesh_fineness_sweep(blade_count=3):
    """
    Mesh independence study: sweep over three mesh fineness settings
    (coarse → medium → fine) for a given blade count.

    Returns dict: label -> (dofs, hemolysis, stagnation, torque)
    """
    mesh_settings = [
        (0.50, 2.20, "coarse"),   # ~89k DOF  — confirmed working
        (0.42, 1.75, "medium"),   # ~207k DOF — confirmed working
        (0.40, 1.65, "fine"),     # ~280-320k DOF (h_min=0.35 gave 424k → OOM)
    ]

    results = {}  # label -> (dofs, hemolysis, stagnation, torque)

    for mesh_min, mesh_max, label in mesh_settings:
        factor = mesh_max / 1.6
        print(f"\n[MESH INDEP] {label} mesh  "
              f"(h_min={mesh_min}, h_max={mesh_max}, CLF={factor:.2f})...")
        build_fluid_mesh_from_step(blade_count, TEMP_MESH,
                                   mesh_min=mesh_min, mesh_max=mesh_max)
        h, s, t, d = solve_cfd_and_get_metrics(TEMP_MESH, blade_count)
        results[label] = (d, h, s, abs(t))
        if os.path.exists(TEMP_MESH):
            os.remove(TEMP_MESH)

    return results


def plot_mesh_fineness(results, blade_count=3):
    """
    Mesh independence plot: 3 subplots (hemolysis, stagnation, torque)
    with DOF count on the X-axis. Shows whether the solution has converged.
    """
    labels = [l for l in ["coarse", "medium", "fine"] if l in results]
    dofs  = [results[l][0] for l in labels]
    hem   = [results[l][1] for l in labels]
    stag  = [results[l][2] for l in labels]
    torq  = [results[l][3] for l in labels]

    # X-tick labels: e.g. "coarse\n12 345 DOF"
    xtick_labels = [
        f"{lbl}\n{d:,} DOF"
        for lbl, d in zip(labels, dofs)
    ]

    # Relative change to finest mesh (convergence check)
    def rel_change(values):
        ref = values[-1] if values[-1] != 0 else 1.0
        return [(v - ref) / ref * 100 for v in values]

    fig, axes = plt.subplots(3, 1, figsize=(9, 14), sharex=True)
    fig.suptitle(
        f'Mesh Independence Study  (N={blade_count} blades)',
        fontweight='bold', fontsize=14, y=0.98
    )

    datasets = [
        (hem,  'Hemolysis: Dissipation Rate (W)', 'tab:red',   'o-'),
        (stag, 'Stagnation: Stagnant Volume (%)',  'tab:blue',  's-'),
        (torq, 'Torque Magnitude (N·mm)',           'tab:green', '^-'),
    ]

    for ax, (vals, ylabel, color, marker) in zip(axes, datasets):
        ax.plot(xtick_labels, vals, marker, color=color,
                linewidth=2.5, markersize=10)
        ax.set_ylabel(ylabel, fontweight='bold', fontsize=11)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.tick_params(axis='x', labelsize=10)

        # Annotate each point with its value
        for xi, (xt, v) in enumerate(zip(xtick_labels, vals)):
            ax.annotate(f'{v:.3g}',
                        xy=(xi, v),
                        xytext=(6, 6), textcoords='offset points',
                        fontsize=9, color=color, fontweight='bold')

        # Right-side secondary axis: % relative to finest mesh
        ax2 = ax.twinx()
        rc = rel_change(vals)
        ax2.plot(xtick_labels, rc, 'D--', color='grey',
                 linewidth=1.2, markersize=6, alpha=0.7,
                 label='Δ vs finest (%)')
        ax2.axhline(0, color='grey', linewidth=0.8, linestyle=':')
        ax2.set_ylabel('Δ vs finest mesh (%)', color='grey', fontsize=9)
        ax2.tick_params(axis='y', labelcolor='grey', labelsize=8)
        ax2.legend(loc='upper right', fontsize=8)

    axes[-1].set_xlabel('Mesh level  (DOF count)', fontweight='bold', fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    p = f"{OUTPUT_DIR}/mesh_independence_N{blade_count}.png"
    plt.savefig(p, dpi=300)
    print(f"\n[SAVED] Mesh independence plot: {p}")


if __name__ == "__main__":
    run_sweep()
    # Mesh independence study for the full rotor (3 blades)
    mesh_results = mesh_fineness_sweep(blade_count=3)
    plot_mesh_fineness(mesh_results, blade_count=3)
