from fenics import *
import numpy as np
import matplotlib.pyplot as plt


tol = 1E-14
kappa = 1
alpha = 1000


mesh = UnitSquareMesh(8, 8)

V = FunctionSpace(mesh, "P", 1)

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)

bx0 = CompiledSubDomain("on_boundary && near(x[0], 0, tol)", tol=tol)
bx1 = CompiledSubDomain("on_boundary && near(x[0], 1, tol)", tol=tol)
by0 = CompiledSubDomain("on_boundary && near(x[1], 0, tol)", tol=tol)
by1 = CompiledSubDomain("on_boundary && near(x[1], 1, tol)", tol=tol)

bx0.mark(boundary_markers, 0)
bx1.mark(boundary_markers, 1)
by0.mark(boundary_markers, 2)
by1.mark(boundary_markers, 3)

u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

bcs = [
    DirichletBC(V, u_D, bx0),
    DirichletBC(V, u_D, bx1),
]

ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)

g = Constant(4)

q = Expression(f"(1 + x[0]*x[0] + 2*x[1]*x[1]) * {alpha}", degree=2)

F = kappa*dot(grad(u), grad(v))*dx + alpha * u * v * ds(2) - f*v*dx - g*v*ds(3) - q*v*ds(2)
a = lhs(F)
L = rhs(F)

# a = kappa*dot(grad(u), grad(v))*dx + alpha * u * v * ds(2)
# L = f*v*dx + g*v*ds(3) + q*v*ds(2)

u = Function(V)
solve(a == L, u, bcs)

plot(u)

error_L2 = errornorm(u_D, u, "L2")

print("error_L2  =", error_L2)

plt.show()
