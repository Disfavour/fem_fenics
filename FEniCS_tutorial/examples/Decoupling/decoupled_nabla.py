from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)

t = 0
T = 1
tau_ = 0.02
time_steps = int(T / tau_)
tau = Constant(tau_)
K = 2

tol = 1e-14
domain_size = 5
mesh_size = 200

a = Constant(1)
gam = Constant(1.4)
alf = 2
bet = 20


mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size, "right")

V = VectorFunctionSpace(mesh, "P", 1)
S = FunctionSpace(mesh, "P", 1)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
subdomain_x.mark(boundaries, 1)
subdomain_y.mark(boundaries, 2)

bcs = [DirichletBC(V.sub(0), Constant(0), boundaries, 1),
       DirichletBC(V.sub(1), Constant(0), boundaries, 2)]

ut = TestFunction(V)
rt = TestFunction(S)
u_ = TrialFunction(V)
r_ = TrialFunction(S)
u = project(Constant((0, 0)), V)
r = project(Expression("1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=1), S)
un = Function(V)
rn = Function(S)
un.assign(u)
rn.assign(r)
uk = Function(V)
rk = Function(S)

F1 = (r_-rn)/tau*rt * dx\
    - dot(r_*uk, grad(rt)) * dx

a1, L1 = lhs(F1), rhs(F1)

F2 = dot((r*u_-rn*un)/tau, ut) * dx\
    - inner(outer(r*uk, u_), nabla_grad(ut)) * dx\
    + dot(grad(a * r**gam), ut) * dx

a2, L2 = lhs(F2), rhs(F2)

Es = []

mass_eq = r * dx
energy_eq = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx
mass = assemble(mass_eq)
energy = assemble(energy_eq)
print(f"conservation: t {t:.3f}\t M {mass:.10f}\t E {energy:.10f}")

Es.append(energy)

for t in [(i+1)*tau_ for i in range(time_steps)]:
    uk.assign(un)
    rk.assign(rn)

    for k in range(K):
        solve(a1 == L1, r)
        solve(a2 == L2, u, bcs)
        uk.assign(u)
        rk.assign(r)

    mass = assemble(mass_eq)
    energy = assemble(energy_eq)
    print(f"conservation: t {t:.3f}\t M {mass:.10f}\t E {energy:.10f}")

    Es.append(energy)

    un.assign(u)
    rn.assign(r)

np.save("decoupled_nabla2.npy", np.array(Es))