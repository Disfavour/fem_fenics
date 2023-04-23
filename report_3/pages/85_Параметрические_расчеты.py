import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Параметрические расчеты
"""

col = st.columns([2, 3])
with col[0]:
    nx1 = st.slider("nx1", 1, 30, 8)
    nx2 = st.slider("nx2", 1, 30, 8)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 2)
    degree = st.slider("degree", 1, 2, 1)
    show_mesh = st.checkbox('визуализация сетки')

T = st.slider("T", 0.01, 10.0, 2.0)
num_steps = st.slider("num_steps", 1, 200, 50)
dt = T / num_steps

step = st.slider("step", 1, num_steps, 1)


# T = 2.0
# num_steps = 50
# dt = T / num_steps

#nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx1, nx2, diagonal)

V = FunctionSpace(mesh, "P", degree)

bc = DirichletBC(V, Constant(0), "on_boundary")

u_0 = Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))", degree=2, a=5)
u_n = interpolate(u_0, V)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

u = Function(V)
t = 0
for n in range(num_steps):
    # Обновляем время
    t += dt
    # Считаем u на следующем шаге
    solve(a == L, u, bc)

    if int(t / dt) == step:
        plot(u)
        if show_mesh:
            plot(mesh)
        with col[1]:
            st.pyplot(plt.gcf())

    # Обновляем результат предыдущего шага
    u_n.assign(u)

