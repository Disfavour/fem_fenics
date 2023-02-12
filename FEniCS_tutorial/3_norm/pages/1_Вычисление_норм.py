import streamlit as st


r"""
# Вычисление норм

### norm(v, norm_type='L2', mesh=None)

Вычисление нормы вектора или функции.

#### Параметры

- v - Vector или Function.
- norm_type - тип нормы.

  Для Vector:
  - "l1" - $L^1 = \displaystyle\sum_i |x_i|$
  - "l2" - $L^2 = \sqrt {\displaystyle\sum_i {|x_i|}^2}$
  - "linf" - $L^\infty = \max |x_i|$
  
  Для Function:
  - "L2" - $L^2 = \sqrt {\displaystyle \int \limits_\Omega u^2 dx}$
  - "H1" - $H^1 = \displaystyle \sqrt {\int \limits_\Omega u^2 + (\nabla u)^2 dx}$
  - "H10" - $H^1_0 = \displaystyle \sqrt {\int \limits_\Omega u^2 + (\nabla u)^2 dx},u = 0$ на $\partial \Omega$

- mesh - сетка для вычисления нормы.
"""

with st.expander("Пример"):
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

        vec = u.vector()

        arr = vec.get_local()

    "### Нормы векторов"

    col = st.columns(2)

    with col[0]:
        "#### FEniCS"
        with st.echo():
            l1 = norm(vec, "l1")
        l1
        with st.echo():
            l2 = norm(vec, "l2")
        l2
        with st.echo():
            linf = norm(vec, "linf")
        linf

    with col[1]:
        "#### Пример реализации"
        with st.echo():
            l1_ = arr.sum()
        l1_
        with st.echo():
            l2_ = (arr ** 2).sum() ** 0.5
        l2_
        with st.echo():
            linf_ = arr.max()
        linf_

    if np.allclose(l1, l1_) and np.allclose(l2, l2_) and np.allclose(linf, linf_):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")

    "### Нормы функций"

    "### $L^2$"
    col = st.columns(2)

    with col[0]:
        with st.echo():
            L2 = norm(u, "L2")
        st.write(L2)

    with col[1]:
        with st.echo():
            from scipy import integrate

            fun = lambda y, x: u(x, y) ** 2
            L2_ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
        st.write(L2_)

    if np.allclose(L2, L2_):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")

    "### $H^1$"
    col = st.columns(2)

    with col[0]:
        with st.echo():
            H_1 = norm(u, "H1")
        st.write(H_1)

    with col[1]:
        with st.echo():
            from scipy import integrate

            h = 1e-6

            fun = lambda y, x: \
                u(x, y) ** 2 \
                + ((u(x + h, y) - u(x - h, y)) / (2 * h)) ** 2 \
                + ((u(x, y + h) - u(x, y - h)) / (2 * h)) ** 2

            H_1_ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
        st.write(H_1_)

    if np.allclose(L2, L2_):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")

    "### $H_0^1$"
    col = st.columns(2)

    with col[0]:
        with st.echo():
            H_1_0 = norm(u, "H10")
        st.write(H_1_0)

    with col[1]:
        with st.echo():
            from scipy import integrate

            h = 1e-6

            u_ = lambda x, y: u(x, y) if h < x < 1 - h and h < y < 1 - h else 0

            fun = lambda y, x:\
                u_(x, y) ** 2\
                + ((u_(x + h, y) - u_(x - h, y)) / (2 * h)) ** 2\
                + ((u_(x, y + h) - u_(x, y - h)) / (2 * h)) ** 2

            H_1_0_ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
        st.write(H_1_0_)

    if np.allclose(H_1_0, H_1_0_):
        st.success("Значения достаточно близки")
    else:
        st.error("Значения различаются")
