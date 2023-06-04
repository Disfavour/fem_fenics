from fenics import *
import matplotlib.pyplot as plt
from ufl import nabla_div


L = 1
W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 10, 3, 3)
V = VectorFunctionSpace(mesh, "P", 1)

bc = DirichletBC(V, Constant((0, 0, 0)), CompiledSubDomain("x[0] < 1e-14"))

def epsilon(u):
    #return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*nabla_grad(u)#epsilon(u)

u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

u = Function(V)
solve(a == L, u, bc)
plot(u)
plt.show()

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d) # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, "P", 1)
von_Mises = project(von_Mises, V)
plot(von_Mises)

# Compute magnitude of displacement
# u_magnitude = sqrt(dot(u, u))
# u_magnitude = project(u_magnitude, V)
# plot(u_magnitude, "Displacement magnitude")
# print("min/max u:",
# u_magnitude.vector().array().min(),
# u_magnitude.vector().array().max())
# plt.show()