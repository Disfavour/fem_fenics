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
  - "H1" - $H^1 = \displaystyle \sqrt {\int \limits_\Omega u^2 + (\nabla u)^2 dx}$
  - "H10" - $H^1_0 = \displaystyle \sqrt {\int \limits_\Omega u^2 + (\nabla u)^2 dx},u = 0$ на $\partial \Omega$
- degree_rise - число на которое увеличится степень интерполированного выражения в пространстве конечных элементов
  относительно uh.
- mesh - сетка для вычисления нормы.

Сначала происходит интерполирование u и uh в общее функциональное пространство, затем вычитание.
"""

with st.expander("Пример"):
    with st.echo():
        from fenics import *
        import numpy as np

        mesh = RectangleMesh(Point(0, 0), Point(1, 1), 1, 1, "left")

        V = FunctionSpace(mesh, "CG", 1)

        g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2, domain=mesh)

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

    "#### $L^2$"
    col = st.columns(3)

    with col[0]:
        with st.echo():
            L2 = errornorm(g, u, "L2")
        st.write(L2)

    with col[1]:
        with st.echo():
            L2_ = assemble((g - u) ** 2 * dx) ** 0.5
        st.write(L2_)

    with col[2]:
        with st.echo():
            from scipy import integrate

            fun = lambda y, x: (u(x, y) - g(x, y)) ** 2
            L2__ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
        st.write(L2__)

    if np.allclose(L2, L2_) and np.allclose(L2, L2__):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")

    "#### $H^1$"
    col = st.columns(3)

    with col[0]:
        with st.echo():
            H_1 = errornorm(g, u, "H1")
        st.write(H_1)

    with col[1]:
        with st.echo():
            H_1_ = assemble(((g - u) ** 2 + grad(g - u) ** 2) * dx) ** 0.5
        st.write(H_1_)

    with col[2]:
        with st.echo():
            from scipy import integrate

            h = 1e-6

            fun = lambda y, x: \
                (g(x, y) - u(x, y)) ** 2 \
                + ((g(x + h, y) - g(x - h, y)) / (2 * h) - (u(x + h, y) - u(x - h, y)) / (2 * h)) ** 2 \
                + ((g(x, y + h) - g(x, y - h)) / (2 * h) - (u(x, y + h) - u(x, y - h)) / (2 * h)) ** 2

            H_1__ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
        st.write(H_1__)

    if np.allclose(H_1, H_1_) and np.allclose(H_1, H_1__):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")
