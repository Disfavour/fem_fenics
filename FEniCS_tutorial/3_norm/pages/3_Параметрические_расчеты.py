import streamlit as st
import matplotlib.pyplot as plt
from fenics import *
import numpy as np


r"""
# Параметрические расчеты

Рассматривается ошибка между решениями краевой задачи Дирихле для уравнения Пуассона на 2-х прямоугольных сетках с
разной степенью конечных элементов.
"""

u1, u2 = None, None

g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

col = st.columns(2)

with col[0]:
    "## 1-й набор параметров `u1`"
    p1 = Point(st.slider("x1", -10.0, 10.0, 0.0, key=301), st.slider("x2", -10.0, 10.0, 0.0, key=302))
    p2 = Point(st.slider("x1", -10.0, 10.0, 1.0, key=303), st.slider("x2", -10.0, 10.0, 1.0, key=304))
    nx1 = st.slider("nx1", 1, 100, 10, key=305)
    nx2 = st.slider("nx2", 1, 100, 10, key=306)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 4, key=315)

    degree = st.slider("degree", 1, 5, 2, key=313)

    mesh = RectangleMesh(p1, p2, nx1, nx2, diagonal)

    plot(mesh)
    st.pyplot(plt.gcf())
    plt.clf()

    f"Количество ячеек ${mesh.num_cells()}$"

    V = FunctionSpace(mesh, "CG", degree)

    bc = DirichletBC(V, g, boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u1 = Function(V)
    solve(a == L, u1, bc)

with col[1]:
    "## 2-й набор параметров `u2`"
    p1 = Point(st.slider("x1", -10.0, 10.0, 0.0, key=307), st.slider("x2", -10.0, 10.0, 0.0, key=308))
    p2 = Point(st.slider("x1", -10.0, 10.0, 1.0, key=309), st.slider("x2", -10.0, 10.0, 1.0, key=310))
    nx1 = st.slider("nx1", 1, 100, 5, key=311)
    nx2 = st.slider("nx2", 1, 100, 5, key=312)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), key=316)

    degree = st.slider("degree", 1, 5, 1, key=314)

    mesh = RectangleMesh(p1, p2, nx1, nx2, diagonal)

    plot(mesh)
    st.pyplot(plt.gcf())
    plt.clf()

    f"Количество ячеек ${mesh.num_cells()}$"

    V = FunctionSpace(mesh, "CG", degree)

    bc = DirichletBC(V, g, boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u2 = Function(V)
    solve(a == L, u2, bc)

f"$L^2 = {errornorm(u1, u2, 'L2')}$ (`errornorm(u1, u2, 'L2')`)"
f"$H^1 = {errornorm(u1, u2, 'H1')}$ (`errornorm(u1, u2, 'H1')`)"
f"$H_0^1 = {errornorm(u1, u2, 'H10')}$ (`errornorm(u1, u2, 'H10')`)"
