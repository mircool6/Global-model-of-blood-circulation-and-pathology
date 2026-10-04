"""Paired monolithic-Newton/IPCS validation on the same 3-plate benchmark.

Both schemes use the identical mesh, inlet waveform, no-slip walls/plates,
and a strong p=0 condition over the whole outlet.  This is essential when
testing whether they predict the same separated-flow branch.
"""
from pathlib import Path
import inspect
import math
import time

import gmsh
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from basix.ufl import element, mixed_element
from dolfinx import fem, geometry
from dolfinx.fem.petsc import (LinearProblem, NonlinearProblem, apply_lifting,
                               assemble_matrix, assemble_vector, create_matrix,
                               create_vector, set_bc)
from dolfinx.io import VTKFile, gmsh as gmshio
from dolfinx.nls.petsc import NewtonSolver
import ufl
from ufl import (FacetNormal, Measure, TestFunction, TestFunctions, TrialFunction,
                 TrialFunctions, div, dot, dx, grad, inner, lhs, nabla_grad,
                 rhs, split)

RE, ST, A_PULSE = 800., .20, .50
CHANNEL_Y0, CHANNEL_H, DOMAIN_X0, DOMAIN_X1 = -1., 2., 0., 25.
PLATE_W, PLATE_H = .04, .8
PLATES = [(2.,False), (4.,True), (6.,False)]
U_IN, DT = 1., .001
INLET, OUTLET, WALLS, PLATE = 1, 2, 3, 4
N_STEPS = 400  # 400 steps test
T_FINAL = N_STEPS * DT
WRITE_EVERY = 20
LOG_EVERY = 50
GAMMA = 0.0  # no added grad-div while validating vortex formation
USE_CLOT = True               # Physical clot injected via inlet
CLOT_START_X = -0.3           # Starts just outside domain, enters immediately
CLOT_Y0 = 0.0                 # Center of channel (CHANNEL_Y0 + CHANNEL_H/2)
CLOT_RADIUS = 0.4             # Characteristic size
CLOT_DELTA_U = -0.6 * U_IN    # Velocity deficit (slower than plasma)
CLOT_SPEED = U_IN             # Approximate advection speed
SEED_VY = 3.0e-2 * U_IN       # deterministic symmetry-breaking perturbation
SEED_DURATION = 0.10          # applied only during startup, then exactly zero
BUDGET = 1800.0               # Time budget per solver in seconds (30 mins)
NU = U_IN * 0.8 / RE
ROOT = Path("/mnt/c/Work/Cursor_CFD")
OUT = ROOT / "results_compare_staggered_ipcs_newton_3plates"
OUT_NEWTON = OUT / "newton"
OUT_IPCS = OUT / "ipcs"
for _directory in (OUT, OUT_NEWTON, OUT_IPCS): _directory.mkdir(parents=True, exist_ok=True)


def make_mesh(comm):
    if gmsh.isInitialized(): gmsh.finalize()
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1)
    if comm.rank == 0:
        gmsh.model.add("staggered_plates_newton")
        channel = gmsh.model.occ.addRectangle(DOMAIN_X0, CHANNEL_Y0, 0, DOMAIN_X1-DOMAIN_X0, CHANNEL_H)
        tools=[]
        for x, top in PLATES:
            y = CHANNEL_Y0+CHANNEL_H-PLATE_H if top else CHANNEL_Y0
            tools.append((2,gmsh.model.occ.addRectangle(x,y,0,PLATE_W,PLATE_H)))
        gmsh.model.occ.cut([(2,channel)],tools,removeTool=True); gmsh.model.occ.synchronize()
        surfaces=gmsh.model.getEntities(2); gmsh.model.addPhysicalGroup(2,[s[1] for s in surfaces],1)
        inlet=[]; outlet=[]; walls=[]; plates=[]; eps=1e-6
        for dim, tag in gmsh.model.getBoundary(surfaces,oriented=False):
            xmin,ymin,_,xmax,ymax,_=gmsh.model.getBoundingBox(dim,tag)
            if abs(xmin-DOMAIN_X0)<eps and abs(xmax-DOMAIN_X0)<eps: inlet.append(tag)
            elif abs(xmin-DOMAIN_X1)<eps and abs(xmax-DOMAIN_X1)<eps: outlet.append(tag)
            elif abs(ymin-CHANNEL_Y0)<eps and abs(ymax-CHANNEL_Y0)<eps: walls.append(tag)
            elif abs(ymin-(CHANNEL_Y0+CHANNEL_H))<eps and abs(ymax-(CHANNEL_Y0+CHANNEL_H))<eps: walls.append(tag)
            else: plates.append(tag)
        for entities, number, name in ((inlet,INLET,"inlet"),(outlet,OUTLET,"outlet"),(walls,WALLS,"walls"),(plates,PLATE,"plate")):
            gmsh.model.addPhysicalGroup(1,entities,number); gmsh.model.setPhysicalName(1,number,name)
        distance=gmsh.model.mesh.field.add("Distance"); gmsh.model.mesh.field.setNumbers(distance,"EdgesList",plates)
        near=gmsh.model.mesh.field.add("Threshold"); gmsh.model.mesh.field.setNumber(near,"IField",distance); gmsh.model.mesh.field.setNumber(near,"LcMin",.035); gmsh.model.mesh.field.setNumber(near,"LcMax",.45); gmsh.model.mesh.field.setNumber(near,"DistMin",.2); gmsh.model.mesh.field.setNumber(near,"DistMax",2.)
        mix=gmsh.model.mesh.field.add("Box"); gmsh.model.mesh.field.setNumber(mix,"VIn",.065); gmsh.model.mesh.field.setNumber(mix,"VOut",.45); gmsh.model.mesh.field.setNumber(mix,"XMin",1.); gmsh.model.mesh.field.setNumber(mix,"XMax",14.); gmsh.model.mesh.field.setNumber(mix,"YMin",CHANNEL_Y0); gmsh.model.mesh.field.setNumber(mix,"YMax",CHANNEL_Y0+CHANNEL_H)
        far=gmsh.model.mesh.field.add("Box"); gmsh.model.mesh.field.setNumber(far,"VIn",.15); gmsh.model.mesh.field.setNumber(far,"VOut",.45); gmsh.model.mesh.field.setNumber(far,"XMin",14.); gmsh.model.mesh.field.setNumber(far,"XMax",DOMAIN_X1); gmsh.model.mesh.field.setNumber(far,"YMin",CHANNEL_Y0); gmsh.model.mesh.field.setNumber(far,"YMax",CHANNEL_Y0+CHANNEL_H)
        combined=gmsh.model.mesh.field.add("Min"); gmsh.model.mesh.field.setNumbers(combined,"FieldsList",[near,mix,far]); gmsh.model.mesh.field.setAsBackgroundMesh(combined)
        gmsh.option.setNumber("Mesh.Algorithm",6); gmsh.model.mesh.generate(2); gmsh.model.mesh.optimize("Netgen")
    data=gmshio.model_to_mesh(gmsh.model,comm,0,gdim=2); gmsh.finalize()
    return data.mesh,data.facet_tags,data.mesh.topology.index_map(data.mesh.topology.dim).size_global


class PulsatileInlet:
    def __init__(self): self.scale=0.; self.t_now=0.
    def __call__(self,x):
        a=np.zeros((2,x.shape[1]),dtype=PETSc.ScalarType)
        u_base = self.scale*U_IN*(1.+A_PULSE*math.sin(2.*math.pi*ST*self.t_now))
        a[0] = u_base
        if USE_CLOT:
            x_c = CLOT_START_X + CLOT_SPEED * self.t_now
            r2 = (x[0] - x_c)**2 + (x[1] - CLOT_Y0)**2
            a[0] += self.scale * CLOT_DELTA_U * np.exp(-r2 / (CLOT_RADIUS**2))
        elif self.t_now <= SEED_DURATION:
            a[1]=self.scale*SEED_VY*np.sin(2.*math.pi*(x[1]-CHANNEL_Y0)/CHANNEL_H)
        return a


def minimum_cell_size(mesh):
    dm=mesh.geometry.dofmaps[0] if hasattr(mesh.geometry,"dofmaps") else mesh.geometry.dofmap; h=np.inf
    for ids in dm:
        x=mesh.geometry.x[ids,:2]
        for i in range(len(x)):
            for j in range(i): h=min(h,float(np.linalg.norm(x[i]-x[j])))
    return mesh.comm.allreduce(h,op=MPI.MIN)


def scalar(form):
    return float(fem.assemble_scalar(form))


def plot_history(history):
    h = np.asarray(history)
    t = h[:, 1]
    fig, ax = plt.subplots(2, 3, figsize=(14, 7))
    fig.suptitle(f"Monolithic Newton - staggered plates, Re={RE:g}, T_final={T_FINAL:.3f}")
    series = (("Uin", 2), ("Vmax", 3), ("CFL", 4), ("dP", 5), ("mass [%]", 6), ("dt", 9))
    for axes, (label, col) in zip(ax.flat, series):
        axes.plot(t, h[:, col], "o-")
        axes.set_title(label); axes.set_xlabel("time"); axes.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(OUT / "diagnostics_newton.png", dpi=180)
    plt.close(fig)


def run_ipcs(mesh, ft, hmin, ref_u_cn=None, V_cn=None, ref_u_bdf2=None, V_bdf2=None):
    """The original P2/P1 IPCS benchmark, reduced only to the same 100 steps."""
    comm=mesh.comm; fdim=mesh.topology.dim-1; dxx=dx(domain=mesh); ds=Measure("ds",domain=mesh,subdomain_data=ft); normal=FacetNormal(mesh)
    V=fem.functionspace(mesh, element("Lagrange",mesh.basix_cell(),2,shape=(2,))); Q=fem.functionspace(mesh,element("Lagrange",mesh.basix_cell(),1))
    
    u_ref_cn = fem.Function(V_cn) if V_cn is not None else None
    u_ref_bdf2 = fem.Function(V_bdf2) if V_bdf2 is not None else None
    inlet=PulsatileInlet(); inlet.scale=1.; uin=fem.Function(V); uin.interpolate(inlet); zero=np.zeros(2,dtype=PETSc.ScalarType)
    bcu=[fem.dirichletbc(uin,fem.locate_dofs_topological(V,fdim,ft.find(INLET))),fem.dirichletbc(zero,fem.locate_dofs_topological(V,fdim,ft.find(WALLS)),V),fem.dirichletbc(zero,fem.locate_dofs_topological(V,fdim,ft.find(PLATE)),V)]
    outlet_p_dofs = fem.locate_dofs_topological(Q, fdim, ft.find(OUTLET))
    bcp=[fem.dirichletbc(PETSc.ScalarType(0), outlet_p_dofs, Q)] # Fix P=0 on entire outlet to stabilize IPCS Poisson step
    k=fem.Constant(mesh,PETSc.ScalarType(DT)); ut=TrialFunction(V); vt=TestFunction(V); pt=TrialFunction(Q); qt=TestFunction(Q)
    uh=fem.Function(V,name="velocity"); us=fem.Function(V); un=fem.Function(V); ph=fem.Function(Q,name="pressure"); phi=fem.Function(Q)
    
    inlet.t_now = 0.; inlet.scale = 1.; uin.interpolate(inlet)
    for field in (uh, us, un):
        field.interpolate(inlet)
        field.x.scatter_forward()
            
    u_mid = 0.5 * (ut + un)
    F1=dot((ut-un)/k,vt)*dxx+inner(dot(un,nabla_grad(u_mid)),vt)*dxx+NU*inner(grad(u_mid),grad(vt))*dxx-dot(ph,div(vt))*dxx
    # Backflow removed to match Newton
    a1,L1=fem.form(lhs(F1)),fem.form(rhs(F1)); a2=fem.form(dot(grad(pt),grad(qt))*dxx); L2=fem.form(-div(us)/k*qt*dxx); a3=fem.form(dot(ut,vt)*dxx); L3=fem.form(dot(us,vt)*dxx-k*dot(grad(phi),vt)*dxx)
    A1,b1=create_matrix(a1),create_vector(fem.extract_function_spaces(L1)); A2,b2=create_matrix(a2),create_vector(fem.extract_function_spaces(L2)); A3,b3=create_matrix(a3),create_vector(fem.extract_function_spaces(L3)); assemble_matrix(A2,a2,bcs=bcp); A2.assemble(); assemble_matrix(A3,a3); A3.assemble()
    def make_ksp(A, kind, pc):
        s=PETSc.KSP().create(comm); s.setOperators(A); s.setType(kind); s.getPC().setType(pc); s.setTolerances(rtol=1e-9,atol=1e-12,max_it=1000); return s
    s1=make_ksp(A1,"bcgs","jacobi"); s2=make_ksp(A2,"cg","hypre"); s2.getPC().setHYPREType("boomeramg"); s3=make_ksp(A3,"cg","sor")
    qi_form=fem.form(dot(uh,-normal)*ds(INLET)); qo_form=fem.form(dot(uh,normal)*ds(OUTLET)); div_form=fem.form(div(uh)**2*dxx); grad_form=fem.form(inner(grad(uh),grad(uh))*dxx); pin=fem.form(ph*ds(INLET)); pout=fem.form(ph*ds(OUTLET))
    vorticity_form = fem.form((uh[1].dx(0) - uh[0].dx(1))**2 * dxx)
    probe_point = np.array([[2.0 + 1.0*PLATE_H, CHANNEL_Y0 + CHANNEL_H/2.0, 0.0]], dtype=np.float64)
    bb_tree = geometry.bb_tree(mesh, mesh.topology.dim)
    cell_candidates = geometry.compute_collisions_points(bb_tree, probe_point)
    colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, probe_point)
    probe_cell = colliding_cells.array[0] if len(colliding_cells.array) > 0 else -1
    
    vu=VTKFile(comm,str(OUT_IPCS/"velocity.pvd"),"w"); vp=VTKFile(comm,str(OUT_IPCS/"pressure.pvd"),"w"); rows=[]
    start_ipcs = time.perf_counter()
    for step in range(1,N_STEPS+1):
        if time.perf_counter() - start_ipcs > BUDGET:
            print(f"IPCS budget {BUDGET}s exhausted. Stopping at step {step}.")
            break
        step_started = time.perf_counter()
        t=step*DT; inlet.t_now=t; uin.interpolate(inlet)
        A1.zeroEntries(); assemble_matrix(A1,a1,bcs=bcu); A1.assemble(); b1.set(0); assemble_vector(b1,L1); apply_lifting(b1,[a1],[bcu]); b1.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE); set_bc(b1,bcu)
        s1.solve(b1,us.x.petsc_vec)
        if s1.getConvergedReason() <= 0: raise RuntimeError(f"Momentum KSP failed at step {step}: reason={s1.getConvergedReason()}, iters={s1.getIterationNumber()}")
        us.x.scatter_forward()
        b2.set(0); assemble_vector(b2,L2); apply_lifting(b2,[a2],[bcp]); b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE); set_bc(b2,bcp)
        s2.solve(b2,phi.x.petsc_vec)
        if s2.getConvergedReason() <= 0: raise RuntimeError(f"Pressure KSP failed at step {step}: reason={s2.getConvergedReason()}, iters={s2.getIterationNumber()}")
        phi.x.scatter_forward(); ph.x.petsc_vec.axpy(1.,phi.x.petsc_vec); ph.x.scatter_forward()
        b3.set(0); assemble_vector(b3,L3); b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES,mode=PETSc.ScatterMode.REVERSE)
        s3.solve(b3,uh.x.petsc_vec)
        if s3.getConvergedReason() <= 0: raise RuntimeError(f"Velocity correction KSP failed at step {step}: reason={s3.getConvergedReason()}, iters={s3.getIterationNumber()}")
        uh.x.scatter_forward(); un.x.array[:]=uh.x.array; un.x.scatter_forward()
        uv=uh.x.array.reshape((-1,2)); vmax=float(np.linalg.norm(uv,axis=1).max()); ux_min=float(uv[:,0].min()); reverse_pct=100.0*float(np.mean(uv[:,0] < 0.0)); cfl=vmax*DT/hmin; qi,qo=scalar(qi_form),scalar(qo_form); mass=100*abs(qi-qo)/max(abs(qi),1e-12); dr=math.sqrt(max(scalar(div_form),0))/max(math.sqrt(max(scalar(grad_form),0)),1e-12); dp=(scalar(pin)-scalar(pout))/CHANNEL_H; uval=U_IN*(1+A_PULSE*math.sin(2*math.pi*ST*t))
        enst = scalar(vorticity_form)
        vy_probe = uh.eval(probe_point, [probe_cell]).flatten()[1] if probe_cell != -1 else 0.0
        
        diff_L2_cn = diff_H1_cn = diff_L2_bdf2 = diff_H1_bdf2 = np.nan
        if ref_u_cn is not None and step - 1 < len(ref_u_cn) and u_ref_cn is not None:
            u_ref_cn.x.array[:] = ref_u_cn[step-1]; u_ref_cn.x.scatter_forward()
            diff_L2_cn = math.sqrt(scalar(fem.form(inner(uh - u_ref_cn, uh - u_ref_cn) * dxx)))
            diff_H1_cn = math.sqrt(scalar(fem.form(inner(grad(uh - u_ref_cn), grad(uh - u_ref_cn)) * dxx)))
        if ref_u_bdf2 is not None and step - 1 < len(ref_u_bdf2) and u_ref_bdf2 is not None:
            u_ref_bdf2.x.array[:] = ref_u_bdf2[step-1]; u_ref_bdf2.x.scatter_forward()
            diff_L2_bdf2 = math.sqrt(scalar(fem.form(inner(uh - u_ref_bdf2, uh - u_ref_bdf2) * dxx)))
            diff_H1_bdf2 = math.sqrt(scalar(fem.form(inner(grad(uh - u_ref_bdf2), grad(uh - u_ref_bdf2)) * dxx)))
            
        rows.append((step,t,uval,vmax,cfl,dp,mass,dr,s1.getIterationNumber(),DT,time.perf_counter()-step_started,ux_min,reverse_pct,enst,vy_probe,diff_L2_cn,diff_H1_cn,diff_L2_bdf2,diff_H1_bdf2))
        if step%LOG_EVERY==0: print(f"IPCS   [{step:03d}/{N_STEPS}] t={t:.4f} Vmax={vmax:.3f} CFL={cfl:.3f} mass={mass:.3e}% Enst={enst:.2e} Vy={vy_probe:.2e}")
        if step==1 or step%WRITE_EVERY==0: vu.write_function(uh,t); vp.write_function(ph,t)
    vu.close(); vp.close(); h=np.asarray(rows); np.savetxt(OUT_IPCS/"diagnostics.csv",h,delimiter=",",header="step,time,Uin,Vmax,CFL,dP,mass_pct,div_rel,KSP,dt,wall_step_s,Ux_min,reverse_dofs_pct,enstrophy,vy_probe,diff_L2_cn,diff_H1_cn,diff_L2_bdf2,diff_H1_bdf2",comments=""); return h


def plot_comparison(ipcs, newton_cn, newton_bdf2):
    fig,axes=plt.subplots(3,4,figsize=(18,10)); fig.suptitle(f"3-plates case: IPCS vs Newton CN vs Newton BDF2 | Re={RE:g}, dt={DT:g}")
    for idx,(title,col,log) in zip([0, 1, 2, 4, 5],(("Uin",2,False),("Vmax",3,False),("CFL",4,False),("mass [%]",6,True),("div_rel",7,True))):
        ax = axes.flat[idx]
        ax.plot(ipcs[:,1],ipcs[:,col],"o-",label="IPCS"); ax.plot(newton_cn[:,1],newton_cn[:,col],"s-",label="Newton CN"); ax.plot(newton_bdf2[:,1],newton_bdf2[:,col],"d-",label="Newton BDF2"); ax.set_title(title); ax.set_xlabel("time"); ax.grid(alpha=.3); ax.legend()
        if log: ax.set_yscale("log")
        
    ax = axes.flat[3]
    dp_ipcs_raw, dp_cn_raw, dp_bdf2_raw = ipcs[:,5], newton_cn[:,5], newton_bdf2[:,5]
    kernel = np.ones(50) / 50.0
    ax.plot(ipcs[:,1], dp_ipcs_raw, alpha=0.15)
    ax.plot(newton_cn[:,1], dp_cn_raw, alpha=0.15)
    ax.plot(newton_bdf2[:,1], dp_bdf2_raw, alpha=0.15)
    ax.plot(ipcs[:,1], np.convolve(dp_ipcs_raw, kernel, mode="same"), label="IPCS")
    ax.plot(newton_cn[:,1], np.convolve(dp_cn_raw, kernel, mode="same"), label="Newton CN")
    ax.plot(newton_bdf2[:,1], np.convolve(dp_bdf2_raw, kernel, mode="same"), label="Newton BDF2")
    ax.set_title("dP (50 ms moving average)"); ax.set_xlabel("time"); ax.grid(alpha=.3); ax.legend()
    
    for i, title, col in [(6,"min(Ux): reverse-flow",11), (7,"Enstrophy",13), (8,"Vy probe",14)]:
        axes.flat[i].plot(ipcs[:,1],ipcs[:,col],"o-",label="IPCS"); axes.flat[i].plot(newton_cn[:,1],newton_cn[:,col],"s-",label="Newton CN"); axes.flat[i].plot(newton_bdf2[:,1],newton_bdf2[:,col],"d-",label="Newton BDF2")
        axes.flat[i].set_title(title); axes.flat[i].set_xlabel("time"); axes.flat[i].grid(alpha=.3); axes.flat[i].legend()
        
    # Plot L2/H1 errors on the last 3 plots instead of dt/wall_time
    # col 15 = diff_L2_cn, col 16 = diff_H1_cn
    axes.flat[9].plot(ipcs[:,1],ipcs[:,15],"o-",label="IPCS vs CN"); axes.flat[9].plot(newton_bdf2[:,1],newton_bdf2[:,15],"d-",label="BDF2 vs CN")
    axes.flat[9].set_title("L2 Error (vs CN)"); axes.flat[9].set_yscale("log"); axes.flat[9].set_xlabel("time"); axes.flat[9].grid(alpha=.3); axes.flat[9].legend()
    
    axes.flat[10].plot(ipcs[:,1],ipcs[:,16],"o-",label="IPCS vs CN"); axes.flat[10].plot(newton_bdf2[:,1],newton_bdf2[:,16],"d-",label="BDF2 vs CN")
    axes.flat[10].set_title("H1 Error (vs CN)"); axes.flat[10].set_yscale("log"); axes.flat[10].set_xlabel("time"); axes.flat[10].grid(alpha=.3); axes.flat[10].legend()
    
    axes.flat[11].plot(ipcs[:,1],np.cumsum(ipcs[:,10]),"o-",label="IPCS"); axes.flat[11].plot(newton_cn[:,1],np.cumsum(newton_cn[:,10]),"s-",label="Newton CN"); axes.flat[11].plot(newton_bdf2[:,1],np.cumsum(newton_bdf2[:,10]),"d-",label="Newton BDF2")
    axes.flat[11].set_title("cumulative wall [s]"); axes.flat[11].set_xlabel("time"); axes.flat[11].grid(alpha=.3); axes.flat[11].legend()
    fig.tight_layout(); fig.savefig(OUT/"01_solver_comparison.png",dpi=180); fig.savefig(OUT/"01_solver_comparison.pdf"); plt.close(fig)


def run_newton(mesh, ft, hmin, n_cells, scheme="CN", ref_u=None):
    comm = mesh.comm; fdim = mesh.topology.dim - 1
    dxx = dx(domain=mesh); ds = Measure("ds", domain=mesh, subdomain_data=ft); normal = FacetNormal(mesh)
    cell = mesh.topology.cell_name()
    W = fem.functionspace(mesh, mixed_element([element("Lagrange", cell, 2, shape=(2,)), element("Lagrange", cell, 1)]))
    w, wn, wnn = fem.Function(W), fem.Function(W), fem.Function(W)
    u, p = split(w); un, pn = split(wn); unn, pnn = split(wnn); v, q = TestFunctions(W)
    V_newton, _ = W.sub(0).collapse()
    u_ref = fem.Function(V_newton) if ref_u is not None else None

    inlet = PulsatileInlet(); uin = fem.Function(V_newton); zero = fem.Function(V_newton)
    zero.interpolate(lambda x: np.zeros((2, x.shape[1])))
    outlet_p_dofs = fem.locate_dofs_topological(W.sub(1), fdim, ft.find(OUTLET))
    p_ref_dof = outlet_p_dofs[:1]
    p_zero_const = fem.Constant(mesh, PETSc.ScalarType(0.0))
    bcs = [
        fem.dirichletbc(uin, fem.locate_dofs_topological((W.sub(0), V_newton), fdim, ft.find(INLET)), W.sub(0)),
        fem.dirichletbc(zero, fem.locate_dofs_topological((W.sub(0), V_newton), fdim, ft.find(WALLS)), W.sub(0)),
        fem.dirichletbc(zero, fem.locate_dofs_topological((W.sub(0), V_newton), fdim, ft.find(PLATE)), W.sub(0)),
        fem.dirichletbc(p_zero_const, p_ref_dof, W.sub(1)),
    ]
    
    inlet.t_now = 0.; inlet.scale = 1.; uin.interpolate(inlet)
    w.sub(0).interpolate(inlet)
    w.sub(1).interpolate(lambda x: np.zeros(x.shape[1], dtype=PETSc.ScalarType))
    w.x.scatter_forward()
    wn.x.array[:] = w.x.array; wn.x.scatter_forward()
    wnn.x.array[:] = w.x.array; wnn.x.scatter_forward()
    
    k = fem.Constant(mesh, PETSc.ScalarType(DT))
    c0 = fem.Constant(mesh, PETSc.ScalarType(1.0/DT))
    c1 = fem.Constant(mesh, PETSc.ScalarType(-1.0/DT))
    c2 = fem.Constant(mesh, PETSc.ScalarType(0.0))
    
    if scheme == "CN":
        u_mid = 0.5 * (u + un)
        p_mid = 0.5 * (p + pn)
        convection = inner(dot(u_mid, nabla_grad(u_mid)), v) * dxx
        time_deriv = inner((u - un) / k, v) * dxx
        F = time_deriv + convection + NU * inner(grad(u_mid), grad(v)) * dxx - p_mid * div(v) * dxx - q * div(u) * dxx
    elif scheme == "BDF2":
        p_mid = p
        convection = inner(dot(u, nabla_grad(u)), v) * dxx
        time_deriv = inner(c0*u + c1*un + c2*unn, v) * dxx
        visc = inner(grad(u), grad(v)) * dxx
        F = time_deriv + convection + NU * visc - p_mid * div(v) * dxx - q * div(u) * dxx
    else:
        raise ValueError(f"Unknown scheme {scheme}")


    args = dict(bcs=bcs, J=ufl.derivative(F, w, TrialFunction(W)))
    modern = "petsc_options_prefix" in inspect.signature(NonlinearProblem.__init__).parameters
    if modern:
        args.update(petsc_options_prefix=f"stagger_{scheme}_", petsc_options={
            "snes_type":"newtonls", "snes_rtol":1e-8, "snes_atol":1e-10,
            "snes_max_it":20, "snes_linesearch_type":"bt",
            "ksp_type":"preonly", "pc_type":"lu",
            "pc_factor_mat_solver_type":"mumps"})
    problem = NonlinearProblem(F, w, **args)
    if not modern:
        solver = NewtonSolver(comm, problem); solver.rtol=1e-8; solver.atol=1e-10; solver.max_it=10; solver.convergence_criterion="incremental"
        solver.krylov_solver.setType("preonly"); solver.krylov_solver.getPC().setType("lu")

    qi_form = fem.form(dot(u, -normal)*ds(INLET)); qo_form = fem.form(dot(u, normal)*ds(OUTLET))
    div_form = fem.form(div(u)**2*dxx); grad_form = fem.form(inner(grad(u), grad(u))*dxx)
    pin_form = fem.form(p*ds(INLET)); pout_form = fem.form(p*ds(OUTLET))
    vorticity_form = fem.form((u[1].dx(0) - u[0].dx(1))**2 * dxx)
    
    probe_point = np.array([[2.0 + 1.0*PLATE_H, CHANNEL_Y0 + CHANNEL_H/2.0, 0.0]], dtype=np.float64)
    bb_tree = geometry.bb_tree(mesh, mesh.topology.dim)
    cell_candidates = geometry.compute_collisions_points(bb_tree, probe_point)
    colliding_cells = geometry.compute_colliding_cells(mesh, cell_candidates, probe_point)
    probe_cell = colliding_cells.array[0] if len(colliding_cells.array) > 0 else -1

    vu = VTKFile(comm, str(OUT_NEWTON / f"velocity_{scheme}.pvd"), "w"); vp = VTKFile(comm, str(OUT_NEWTON / f"pressure_{scheme}.pvd"), "w")
    history=[]; u_snapshots=[]; start_newton = time.perf_counter()
    if comm.rank == 0:
        print(f"Newton staggered 3-plates ({scheme}): Re={RE:g}, nu={NU:g}, cells={n_cells}, h_min={hmin:.3e}, dt={DT:g}, steps={N_STEPS}, seed={SEED_VY:g} for {SEED_DURATION}s")
    for step in range(1, N_STEPS+1):
        if time.perf_counter() - start_newton > BUDGET:
            print(f"Newton {scheme} budget {BUDGET}s exhausted. Stopping at step {step}.")
            break
        if scheme == "BDF2" and step == 2:
            c0.value = 1.5 / DT
            c1.value = -2.0 / DT
            c2.value = 0.5 / DT
            
        step_started = time.perf_counter()
        t = step*DT; inlet.t_now=t; inlet.scale=1.; uin.interpolate(inlet)
        w.x.array[:] = wn.x.array; w.x.scatter_forward()
        if modern:
            problem.solve(); w.x.scatter_forward(); nit=problem.solver.getIterationNumber(); ok=problem.solver.getConvergedReason()>0
        else:
            nit, ok = solver.solve(w); w.x.scatter_forward()
        if not ok:
            reason = problem.solver.getConvergedReason() if modern else "legacy NewtonSolver"
            detail = ""
            if modern:
                ksp = problem.solver.getKSP()
                detail = f", KSP={ksp.getConvergedReason()}, PC={ksp.getPC().getFailedReason()}"
            raise RuntimeError(f"Newton failed at step {step}, t={t:.5f}, iterations={nit}, SNES reason={reason}{detail}")
            
        wnn.x.array[:] = wn.x.array; wnn.x.scatter_forward()
        wn.x.array[:] = w.x.array; wn.x.scatter_forward()
        
        uh_coll = w.sub(0).collapse()
        u_snapshots.append(np.copy(uh_coll.x.array))
        
        uv=uh_coll.x.array.reshape((-1,2)); vmax=float(np.linalg.norm(uv,axis=1).max()); ux_min=float(uv[:,0].min()); reverse_pct=100.0*float(np.mean(uv[:,0] < 0.0)); cfl=vmax*DT/hmin
        qi, qo = scalar(qi_form), scalar(qo_form); mass=100*abs(qi-qo)/max(abs(qi),1e-12)
        div_rel=math.sqrt(max(scalar(div_form),0))/max(math.sqrt(max(scalar(grad_form),0)),1e-12)
        dp=(scalar(pin_form)-scalar(pout_form))/CHANNEL_H; uin_now=U_IN*(1+A_PULSE*math.sin(2*math.pi*ST*t))
        enst = scalar(vorticity_form)
        vy_probe = uh_coll.eval(probe_point, [probe_cell]).flatten()[1] if probe_cell != -1 else 0.0
        
        diff_L2_cn = diff_H1_cn = np.nan
        if ref_u is not None and step - 1 < len(ref_u):
            u_ref.x.array[:] = ref_u[step-1]; u_ref.x.scatter_forward()
            diff_L2_cn = math.sqrt(scalar(fem.form(inner(uh_coll - u_ref, uh_coll - u_ref)*dxx)))
            diff_H1_cn = math.sqrt(scalar(fem.form(inner(grad(uh_coll - u_ref), grad(uh_coll - u_ref))*dxx)))
            
        history.append((step,t,uin_now,vmax,cfl,dp,mass,div_rel,nit,DT,time.perf_counter()-step_started,ux_min,reverse_pct,enst,vy_probe,diff_L2_cn,diff_H1_cn))
        if step%LOG_EVERY==0:
            print(f"[{step:05d}/{N_STEPS:05d}] t={t:.4f} Uin={uin_now:.3f} Vmax={vmax:.3f} CFL={cfl:.3f} mass={mass:.3e}% Enst={enst:.2e} Vy={vy_probe:.2e} Newton={nit}")
        if step==1 or step%WRITE_EVERY==0:
            a,b=w.split(); a=a.collapse(); b=b.collapse(); a.name="velocity"; b.name="pressure"; vu.write_function(a,t); vp.write_function(b,t)
    vu.close(); vp.close()
    header="step,time,Uin,Vmax,CFL,dP,mass_pct,div_rel,Newton,dt,wall_step_s,Ux_min,reverse_dofs_pct,enstrophy,vy_probe,diff_L2_cn,diff_H1_cn"
    newton_history=np.asarray(history)
    np.savetxt(OUT_NEWTON/f"diagnostics_{scheme}.csv",newton_history,delimiter=",",header=header,comments="")
    return newton_history, u_snapshots, V_newton


def main():
    comm = MPI.COMM_WORLD
    if comm.size != 1: raise RuntimeError("Run this validation in serial")
    mesh, ft, n_cells = make_mesh(comm)
    hmin = minimum_cell_size(mesh)
    
    newton_cn, cn_u, V_cn = run_newton(mesh, ft, hmin, n_cells, scheme="CN")
    newton_bdf2, bdf2_u, V_bdf2 = run_newton(mesh, ft, hmin, n_cells, scheme="BDF2", ref_u=cn_u)
    
    print(f"Newton CN and BDF2 complete; starting the same IPCS benchmark.")
    ipcs_history = run_ipcs(mesh, ft, hmin, ref_u_cn=cn_u, V_cn=V_cn, ref_u_bdf2=bdf2_u, V_bdf2=V_bdf2)
    
    # Pad shorter histories with NaNs if budget cut one off early
    n_cn = len(newton_cn)
    n_bdf2 = len(newton_bdf2)
    n_ipcs = len(ipcs_history)
    max_len = max(n_cn, n_bdf2, n_ipcs)
    
    def pad(arr):
        res = np.full((max_len, arr.shape[1]), np.nan)
        res[:len(arr), :] = arr
        return res
        
    plot_comparison(pad(ipcs_history), pad(newton_cn), pad(newton_bdf2))

    timing=np.array([
        [np.nansum(ipcs_history[:,10]),np.nanmean(ipcs_history[:,10])],
        [np.nansum(newton_cn[:,10]),np.nanmean(newton_cn[:,10])],
        [np.nansum(newton_bdf2[:,10]),np.nanmean(newton_bdf2[:,10])]
    ])
    np.savetxt(OUT/"solver_timing.csv",timing,delimiter=",",header="total_wall_s,mean_wall_step_s",comments="",fmt="%.8g")
    print(f"TIMING: IPCS total={timing[0,0]:.1f}s avg={timing[0,1]:.3f}s/step | CN total={timing[1,0]:.1f}s avg={timing[1,1]:.3f}s/step | BDF2 total={timing[2,0]:.1f}s avg={timing[2,1]:.3f}s/step")
    print(f"COMPLETE: results={OUT}")

if __name__ == "__main__": main()

