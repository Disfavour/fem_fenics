import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Параметрические расчеты
"""

cols = st.columns(2)
with cols[0]:
    k_0 = st.slider("k_0", -5.0, 5.0, 1.0)
with cols[1]:
    k_1 = st.slider("k_0", -5.0, 5.0, 0.01)

cols = st.columns(2)
with cols[0]:
    nx1 = st.slider("nx1", 1, 20, 10)
    nx2 = st.slider("nx2", 1, 20, 10)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 4)
    degree = st.slider("degree", 1, 2, 1)
with cols[1]:
    nx1_2 = st.slider("nx1", 20, 50, 50, key=100)
    nx2_2 = st.slider("nx2", 20, 50, 50, key=101)
    diagonal_2 = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 4, key=102)
    degree_2 = st.slider("degree", 1, 2, 1, key=103)
show_mesh = st.checkbox('визуализация сетки')


class K(UserExpression):
    def __init__(self, materials, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.materials = materials
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.materials[cell.index] == 0:
            values[0] = self.k_0
        else:
            values[0] = self.k_1


def get_u(nx1, nx2, diagonal, degree):
    tol = 1E-14

    mesh = UnitSquareMesh(nx1, nx2, diagonal)

    V = FunctionSpace(mesh, "P", degree)

    materials = MeshFunction("size_t", mesh, mesh.topology().dim())

    subdomain_0 = CompiledSubDomain("x[1] <= 0.5 + tol", tol=tol)
    subdomain_1 = CompiledSubDomain("x[1] >= 0.5 - tol", tol=tol)

    materials.set_all(0)
    subdomain_1.mark(materials, 1)

    kappa = K(materials, k_0, k_1, degree=0)

    bcs = DirichletBC(V, Constant(0), "on_boundary")

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = kappa * dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u = Function(V)
    solve(a == L, u, bcs)

    return u, mesh


cols = st.columns(2)
with cols[0]:
    u, mesh = get_u(nx1, nx2, diagonal, degree)
    plt.colorbar(plot(u))
    if show_mesh:
        plot(mesh)
    st.pyplot(plt.gcf())
    plt.clf()
with cols[1]:
    u2, mesh2 = get_u(nx1_2, nx2_2, diagonal_2, degree_2)
    plt.colorbar(plot(u2))
    if show_mesh:
        plot(mesh2)
    st.pyplot(plt.gcf())
    plt.clf()

with st.columns([1, 2, 1])[1]:
    f"""
    $\| u_e - u \|_2 = {errornorm(u2, u)}$
    """
