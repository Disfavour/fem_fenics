import streamlit as st

r"""
# Параметрические расчеты
"""

col = st.columns([2, 3])
with col[0]:
    nx1 = st.slider("nx1", 1, 20, 8)
    nx2 = st.slider("nx2", 1, 20, 8)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 2)
    degree = st.slider("degree", 1, 2, 1)
    show_mesh = st.checkbox('визуализация сетки')


from fenics import *
import matplotlib.pyplot as plt

u_e = Expression("1 + x[0] + 2 * x[1]", degree=2)

mesh = UnitSquareMesh(nx1, nx2, diagonal)
V = FunctionSpace(mesh, 'P', degree)

bc = DirichletBC(V, u_e, "on_boundary")

u = Function(V)
v = TestFunction(V)
f = Expression("-10 - 10*x[0] - 20*x[1]", degree=2)


def q(u):
    return 1 + u ** 2


F = (q(u) * dot(grad(u), grad(v)) - f * v) * dx
#F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx

solve(F == 0, u, bc)
#u_e = interpolate(u_e, V)
error_L2 = errornorm(u_e, u)

plot(u)
if show_mesh:
    plot(mesh)

plt.title("$u$")
plt.show()

with col[1]:
    st.pyplot(plt.gcf())

with st.columns([1, 2, 1])[1]:
    f"""
    $\| u_e - u \|_2 = {error_L2}$
    """
