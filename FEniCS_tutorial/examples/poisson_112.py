from fenics import *
import numpy as np
import matplotlib.pyplot as plt


mesh = UnitSquareMesh(8, 8)

V = FunctionSpace(mesh, "P", 1)

tol = 1E-14

u_L = Expression("1 + 2*x[1]*x[1]", degree=2)


def boundary_L(x, on_boundary):
    return on_boundary and near(x[0], 0, tol)


bc_L = DirichletBC(V, u_L, boundary_L)

u_R = Expression("2 + 2*x[1]*x[1]", degree=2)


def boundary_R(x, on_boundary):
    return on_boundary and near(x[0], 1, tol)


bc_R = DirichletBC(V, u_R, boundary_R)

bcs = [bc_L, bc_R]

# exact solution
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
g = Expression("-4*x[1]", degree=1)
L = f*v*dx - g*v*ds

u = Function(V)
solve(a == L, u, bcs)

plt.colorbar(plot(u))
plot(mesh)

# vtkfile = File("poisson/solution.pvd")
# vtkfile << u

error_L2 = errornorm(u_D, u, "L2")

vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

print("error_L2  =", error_L2)
print("error_max =", error_max)

plt.show()
