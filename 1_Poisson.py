from fenics import *
import numpy as np
import matplotlib.pyplot as plt
#from dolfin import RectangleMesh


# Создание сетки и функционального пространства
#mesh = UnitSquareMesh(8, 8, "right/left")
mesh = RectangleMesh(Point(0, 0), Point(3, 2), 6, 10, "left")
# plot()
# plt.show()
# plt.clf()

V = FunctionSpace(mesh, "P", 1)

# Граничные условия
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Вариационная постановка
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Вычисление
u = Function(V)
solve(a == L, u, bc)

# Визуализация сетки и решения

plt.colorbar(plot(u))
plot(mesh)

# Сохранение решения в VTK формате
# vtkfile = File("poisson/solution.pvd")
# vtkfile << u

# Вычисление ошибки в L2 норме
error_L2 = errornorm(u_D, u, "L2")

# Вычисление максимальной ошибки на вершинах
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

print("error_L2 =", error_L2)
print("error_max =", error_max)

plt.show()
