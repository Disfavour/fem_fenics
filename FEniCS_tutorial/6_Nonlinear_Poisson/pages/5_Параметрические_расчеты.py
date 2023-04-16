import streamlit as st
import sympy
from fenics import *
import matplotlib.pyplot as plt


r"""
# Параметрические расчеты
"""

u_e_str = st.text_input('Решение', '1 + x + 2*y')

#u_e = lambda x, y: eval(u_e_str)
u_e = Expression(u_e_str.replace("x", "x[0]").replace("y", "x[1]"), degree=3)

def q(u):
    return 1 + u**2

x, y = sympy.symbols("x[0], x[1]")

u = eval(u_e_str)
#u = 1 + x + 2*y

f = - sympy.diff(q(u)*sympy.diff(u, x), x) - sympy.diff(q(u)*sympy.diff(u, y), y)
f = sympy.simplify(f)
u_code = sympy.printing.ccode(u)
f_code = sympy.printing.ccode(f)

u_D = Expression(u_code, degree=1)
f = Expression(f_code, degree=1)

col = st.columns([2, 3])
with col[0]:
    nx1 = st.slider("nx1", 1, 20, 8, key=10)
    nx2 = st.slider("nx2", 1, 20, 8, key=11)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 2, key=12)
    degree = st.slider("degree", 1, 2, 1, key=13)
    show_mesh = st.checkbox('визуализация сетки')

mesh = UnitSquareMesh(nx1, nx2, diagonal)
V = FunctionSpace(mesh, 'P', degree)
u_D = Expression(u_code, degree=1)

bc = DirichletBC(V, u_D, "on_boundary")
u = Function(V)
v = TestFunction(V)
f = Expression(f_code, degree=1)
F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx
solve(F == 0, u, bc)

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