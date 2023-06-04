from fenics import *


set_log_level(LogLevel.WARNING)

t = 0
T = 5
time_steps = 10
tau_ = T / time_steps
tau = Constant(tau_)

tol = 1e-14
size = 5
mesh_size = 200

a = Constant(1)
gam = Constant(1.4)
alf = 2
bet = 20


mesh = RectangleMesh(Point(-size, -size), Point(size, size), mesh_size, mesh_size, "right")

V = VectorElement("P", mesh.ufl_cell(), 1)      # triangle
S = FiniteElement("P", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([V, S]))   # V * S

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=size)
subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=size)
subdomain_x.mark(boundaries, 0)
subdomain_y.mark(boundaries, 1)

bcs = [DirichletBC(W.sub(0).sub(0), Constant(0), boundaries, 0),
       DirichletBC(W.sub(0).sub(1), Constant(0), boundaries, 1)]

ut, rt = TestFunctions(W)
w = project(Expression(("0", "0", "1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))"), alf=alf, bet=bet, degree=1), W)
wn = Function(W)
wn.assign(w)
u, r = split(w)
un, rn = split(wn)

F = (r-rn)/tau * rt * dx\
    - dot(r*u, grad(rt)) * dx\
    + dot((r*u-rn*un)/tau, ut) * dx\
    - inner(outer(r*u, u), grad(ut)) * dx\
    + dot(grad(a * r**gam), ut) * dx

sm = r * dx
se = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx
sM = assemble(sm)
sE = assemble(se)
print(f"conservation: t {t:.3f}\t M {sM:.10f}\t E {sE:.10f}")

for t in [(i+1)*tau_ for i in range(time_steps)]:
    J = derivative(F, w)
    solve(F == 0, w, bcs, J=J, solver_parameters={"newton_solver":
                                                 {"relative_tolerance": tol, "absolute_tolerance": tol,
                                                  "relaxation_parameter": 1.0}})

    sM = assemble(sm)
    sE = assemble(se)
    print(f"conservation: t {t:.3f}\t M {sM:.10f}\t E {sE:.10f}")

    wn.assign(w)
