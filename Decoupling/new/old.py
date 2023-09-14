from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)

ts = []
ms = []
Es = []

mesh_size = 100
degree = 1
tau = 0.01
T = 5
t = 0
time_steps = [(i + 1) * tau for i in range(int(T / tau))]

domain_size = 5

alf = 2
bet = 20

tau_ = Constant(tau)
a = Constant(1)
gam = Constant(2)

mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size,
                     "right")

V = VectorElement("P", mesh.ufl_cell(), degree)
S = FiniteElement("P", mesh.ufl_cell(), degree)
W = FunctionSpace(mesh, MixedElement([V, S]))

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
subdomain_x.mark(boundaries, 1)
subdomain_y.mark(boundaries, 2)

bcs = [DirichletBC(W.sub(0).sub(0), Constant(0), boundaries, 1),
       DirichletBC(W.sub(0).sub(1), Constant(0), boundaries, 2)]

ut, rt = TestFunctions(W)
w = Function(W)
wn = Function(W)
u, r = split(w)
un, rn = split(wn)

F = (r-rn)/tau_ * rt * dx \
    - dot(r*u, grad(rt)) * dx \
    + dot((r*u-rn*un)/tau_, ut) * dx \
    - inner(outer(r*u, u), nabla_grad(ut)) * dx \
    + dot(grad(a * r**gam), ut) * dx

m_eq = r * dx
E_eq = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx

# t = 0
w.assign(project(Expression(("0", "0", "1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))"), alf=alf, bet=bet, degree=degree), W))

m, E = assemble(m_eq), assemble(E_eq)
print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

ts.append(t)
ms.append(m)
Es.append(E)

wn.assign(w)

for t in time_steps:
    solve(F == 0, w, bcs)

    m, E = assemble(m_eq), assemble(E_eq)
    print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

    ts.append(t)
    ms.append(m)
    Es.append(E)

    wn.assign(w)

ts = np.array(ts)
ms = np.array(ms)
Es = np.array(Es)

np.save(f"data/old.npy", np.array((ts, ms, Es)).T)
