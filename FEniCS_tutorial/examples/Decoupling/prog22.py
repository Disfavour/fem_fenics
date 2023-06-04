from fenics import *


set_log_level(LogLevel.WARNING)

t = 0
T = 1
time_steps = 35
tau_ = T / time_steps
tau = Constant(tau_)
K = 1

tol = 1e-14
size = 5
mesh_size = 200

a = Constant(1)
gam = Constant(1.4)
alf = 2
bet = 20


mesh = RectangleMesh(Point(-size, -size), Point(size, size), mesh_size, mesh_size, "right")

V = VectorFunctionSpace(mesh, "P", 1)
S = FunctionSpace(mesh, "P", 1)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=size)
subdomain_x.mark(boundaries, 0)
subdomain_y.mark(boundaries, 1)

bcs = [DirichletBC(V.sub(0), Constant(0), boundaries, 0),
       DirichletBC(V.sub(1), Constant(0), boundaries, 1)]

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
    - inner(outer(r*u_, uk), grad(ut)) * dx\
    + dot(grad(a * r**gam), ut) * dx

a2, L2 = lhs(F2), rhs(F2)

sm = r * dx
se = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx
sM = assemble(sm)
sE = assemble(se)
print(f"conservation: t {t:.3f}\t M {sM:.10f}\t E {sE:.10f}")

for t in [(i+1)*tau_ for i in range(time_steps)]:
    uk.assign(un)
    rk.assign(rn)

    for k in range(K):
        solve(a1 == L1, r)
        solve(a2 == L2, u, bcs)
        uk.assign(u)
        rk.assign(r)

    sM = assemble(sm)
    sE = assemble(se)
    print(f"conservation: t {t:.3f}\t M {sM:.10f}\t E {sE:.10f}")

    un.assign(u)
    rn.assign(r)
