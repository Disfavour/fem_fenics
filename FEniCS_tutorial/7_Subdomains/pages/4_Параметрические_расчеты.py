import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


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
k_0 = st.slider("k_0", -5.0, 5.0, 1.0)
k_1 = st.slider("k_0", -5.0, 5.0, 0.01)

from fenics import *
import matplotlib.pyplot as plt

tol = 1E-14
# k_0 = 1.0
# k_1 = 0.01

mesh = UnitSquareMesh(nx1, nx2, diagonal)

V = FunctionSpace(mesh, "P", degree)

materials = MeshFunction("size_t", mesh, mesh.topology().dim())

subdomain_0 = CompiledSubDomain("x[1] <= 0.5 + tol", tol=tol)
subdomain_1 = CompiledSubDomain("x[1] >= 0.5 - tol", tol=tol)

materials.set_all(0)
subdomain_1.mark(materials, 1)


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


kappa = K(materials, k_0, k_1, degree=0)


bcs = DirichletBC(V, Constant(0), "on_boundary")

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = kappa * dot(grad(u), grad(v)) * dx
L = f * v * dx

u = Function(V)
solve(a == L, u, bcs)

plt.colorbar(plot(u))
if show_mesh:
    plot(mesh)

with col[1]:
    st.pyplot(plt.gcf())
