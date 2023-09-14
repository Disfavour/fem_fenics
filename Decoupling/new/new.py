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

W = VectorElement("P", mesh.ufl_cell(), degree)
S = FiniteElement("P", mesh.ufl_cell(), degree)
V = FunctionSpace(mesh, MixedElement([W, S]))

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
subdomain_x.mark(boundaries, 1)
subdomain_y.mark(boundaries, 2)

bcs = [DirichletBC(V.sub(0).sub(0), Constant(0), boundaries, 1),
       DirichletBC(V.sub(0).sub(1), Constant(0), boundaries, 2)]

wt, st = TestFunctions(V)
v = Function(V)
vn = Function(V)
w, s = split(v)
wn, sn = split(vn)


def wm():
    return (w + wn) / 2


def sm():
    return (s + sn) / 2


def um():
    return wm() / sm()


F = (s - sn) / tau_ * st * dx \
    + (div(wm()) + dot(um(), grad(sm()))) / 2 * st * dx \
    + dot((w - wn) / tau_, wt) * dx \
    + dot(div(outer(um(), wm())) + dot(um(), grad(wm())) / 2, wt) * dx \
    + dot(1 / sm() * grad(a * sm() ** (2 * gam)), wt) * dx

m_eq = s**2 * dx
E_eq = (dot(w, w) / 2 + a * s ** (2*gam) / (gam - 1)) * dx

# t = 0
v.assign(project(Expression(("0", "0", "pow(1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])), 0.5)"), alf=alf, bet=bet, degree=degree), V))

t1 = project(Expression(("0", "0", "1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))"), alf=alf, bet=bet, degree=degree), V)
t2 = project(Expression(("0", "0", "pow(1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])), 0.5)"), alf=alf, bet=bet, degree=degree), V)

t11, t12 = split(t1)
t21, t22 = split(t2)

print(assemble((dot(w, w) / 2 + a * t12 ** gam / (gam - 1)) * dx), "ga")
print(assemble((dot(w, w) / 2 + a * ((t22 * t22) ** gam)) / (gam - 1) * dx), "ga")
print(assemble((dot(w, w) / 2 + a * t22 ** (2*gam) / (gam - 1)) * dx), "ga")

m, E = assemble(m_eq), assemble(E_eq)
print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

ts.append(t)
ms.append(m)
Es.append(E)

vn.assign(v)

for t in time_steps:
    solve(F == 0, v, bcs)

    m, E = assemble(m_eq), assemble(E_eq)
    print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

    ts.append(t)
    ms.append(m)
    Es.append(E)

    vn.assign(v)

ts = np.array(ts)
ms = np.array(ms)
Es = np.array(Es)

np.save(f"data/new.npy", np.array((ts, ms, Es)).T)
