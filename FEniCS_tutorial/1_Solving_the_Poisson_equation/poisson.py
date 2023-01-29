import dolfin
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
#from dolfin import RectangleMesh
import gmsh
from scipy import integrate



# Создание сетки и функционального пространства
#mesh = UnitSquareMesh(8, 8, "right/left")
mesh = RectangleMesh(Point(0, 0), Point(1, 1), 1, 1, "left")


#plot(Mesh().coordinates() = np.array([[0], [1]]))
#plt.show()
#plt.clf()

V = FunctionSpace(mesh, "P", 1)

# Граничные условия
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, g, boundary)

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
error_L2 = errornorm(g, u, "L2")

n1 = norm(u)
print(n1)
fun = lambda y, x: u(x, y) ** 2
n2 = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
print(n2)

print(error_L2)
fun = lambda y, x: (g(x, y) - u(x, y)) ** 2
err2 = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
print(err2)

nodal_values_u = u.vector()
array_u = nodal_values_u.get_local()
print(array_u)


print(norm(u.vector(), "l1"))
print(sum(array_u))

print(norm(u.vector(), "l2"))
print(sum([i ** 2 for i in array_u]) ** 0.5)

print(norm(u.vector(), "linf"))
print(max(array_u))

# Вычисление максимальной ошибки на вершинах
vertex_values_g = g.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_g - vertex_values_u))

print("error_L2 =", error_L2)
print("error_max =", error_max)

print(type(u))
print(type(g))
print(type(g - u))

print(errornorm(g, u, "L2"))
print(errornorm(g, u, "L2", mesh=mesh))

print(norm(u, "L2"))

print(type(u.vector()))

print(norm(u.vector(), "linf"))

print(type(g))

print(type(interpolate(g, V) - u))


#plt.show()
