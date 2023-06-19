from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)

t = 0
T = 5
tau_ = 0.005
time_steps = int(T / tau_)
tau = Constant(tau_)
K = 5

tol = 1e-14
domain_size = 5
mesh_size = 200

a = Constant(1)
gam = Constant(1.4)
alf = 2
bet = 20


mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size, "right")

U = FunctionSpace(mesh, "P", 1)
V = FunctionSpace(mesh, "P", 1)
R = FunctionSpace(mesh, "P", 1)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
subdomain_x.mark(boundaries, 1)
subdomain_y.mark(boundaries, 2)

bc_u = DirichletBC(U, Constant(0), boundaries, 1)
bc_v = DirichletBC(V, Constant(0), boundaries, 2)

ut, vt, rt = TestFunction(U), TestFunction(V), TestFunction(R)
u_, v_, r_ = TrialFunction(U), TrialFunction(V), TrialFunction(R)
u, v, r = Function(U), Function(V), Function(R)
un, vn, rn = Function(U), Function(V), Function(R)
uk, vk, rk = Function(U), Function(V), Function(R)

un = project(Constant(0), U)
vn = project(Constant(0), V)
rn = project(Expression("1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=1), R)

u.assign(un)
v.assign(vn)
r.assign(rn)

F1 = (r_-rn)/tau * rt * dx\
    + (r_*uk).dx(0) * rt * dx + (r_*vk).dx(1) * rt * dx

a1, L1 = lhs(F1), rhs(F1)

F2 = (r*u_ - rn*un)/tau * ut * dx\
     + (r*u_*uk).dx(0) * ut * dx + (r*u_*vk).dx(1) * ut * dx\
     + gam * r ** (gam - 1) * r.dx(0) * ut * dx
a2, L2 = lhs(F2), rhs(F2)

F3 = (r*v_ - rn*vn)/tau * vt * dx\
     + (r*v_*u).dx(0) * vt * dx + (r*v_*vk).dx(1) * vt * dx\
     + gam * r ** (gam - 1) * r.dx(1) * vt * dx
a3, L3 = lhs(F3), rhs(F3)

mass_eq = r * dx
energy_eq = 0.5 * r * (u * u + v * v) * dx + 1. / (gam - 1) * r ** gam * dx

mass = assemble(mass_eq)
energy = assemble(energy_eq)
print(f"conservation: t {t:.3f}\t M {mass:.10f}\t E {energy:.10f}")

Es = [energy]

for t in [(i+1)*tau_ for i in range(time_steps)]:
    uk.assign(un)
    vk.assign(vn)
    rk.assign(rn)

    for k in range(K):
        solve(a1 == L1, r)
        solve(a2 == L2, u, bc_u)
        solve(a3 == L3, v, bc_v)

        uk.assign(u)
        vk.assign(v)
        rk.assign(r)

    mass = assemble(mass_eq)
    energy = assemble(energy_eq)
    print(f"conservation: t {t:.3f}\t M {mass:.10f}\t E {energy:.10f}")

    Es.append(energy)

    un.assign(u)
    vn.assign(v)
    rn.assign(r)

np.save("data/seidel5.npy", np.array(Es))
