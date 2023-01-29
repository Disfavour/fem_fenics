import streamlit as st


r"""
# Вычисление ошибок

### errornorm(u, uh, norm_type="l2", degree_rise=3, mesh=None)

Вычисление ошибки `u - uh` в выбранной норме.

#### Параметры

- u - Function.
- uh - Function.
- norm_type - тип нормы.
  - "L2" - $L^2 = \sqrt {\displaystyle \int \limits_\Omega (u - uh)^2 dx}$
  - "H1" - $H^1$
  - "H10" - $H^1_0$
- degree_rise - число на которое увеличится степень интерполированного выражения в пространстве конечных элементов
  относительно uh.
- mesh - сетка для вычисления нормы.

Сначала происходит интерполирование u и uh в общее функциональное пространство, затем вычитание.
"""

with st.echo():
    from fenics import *
    import numpy as np

    mesh = RectangleMesh(Point(0, 0), Point(1, 1), 1, 1, "left")

    V = FunctionSpace(mesh, "CG", 1)

    g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, g, boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u = Function(V)
    solve(a == L, u, bc)

col = st.columns(2)

with col[0]:
    "#### FEniCS"
    with st.echo():
        L2 = errornorm(g, u, "L2")
    L2


with col[1]:
    "#### Пример реализации"
    with st.echo():
        from scipy import integrate

        fun = lambda y, x: (g(x, y) - u(x, y)) ** 2
        L2_ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
    L2_

if np.allclose(L2, L2_):
    st.success("Значения достаточно близки")
else:
    st.error("Значения различаются")
