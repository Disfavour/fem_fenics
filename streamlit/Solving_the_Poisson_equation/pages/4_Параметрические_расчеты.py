import streamlit as st


r"""
# Параметрические расчеты
"""

"$n$"
mesh_size_x_0 = st.slider("", 1, 100, 8, key=1)
"$m$"
mesh_size_x_1 = st.slider("", 1, 100, 8, key=2)
degree = st.slider("Степень конечного элемента", 1, 10, 1)

mesh_vis = st.checkbox("Визуализация сетки")


from fenics import *
import matplotlib.pyplot as plt


# Создание сетки и функционального пространства
mesh = UnitSquareMesh(mesh_size_x_0, mesh_size_x_1)
V = FunctionSpace(mesh, "P", degree)

# Граничные условия
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, g, boundary)

# Вариационная постановка
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(grad(u), grad(v))*dx
l = f*v*dx

# Вычисление
u = Function(V)
solve(a == l, u, bc)

error_L2 = errornorm(g, u, "L2")
f"$L^2$ норма $={error_L2}$"

vertex_values_u_D = g.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
f"Максимальная ошибка $={error_max}$"

plt.colorbar(plot(g - u, title="g - u"))
if mesh_vis:
    plot(mesh)
st.pyplot(plt.gcf())


plt.clf()
# Визуализация сетки и решения
plt.colorbar(plot(u, title="u"))
if mesh_vis:
    plot(mesh)
st.pyplot(plt.gcf())
