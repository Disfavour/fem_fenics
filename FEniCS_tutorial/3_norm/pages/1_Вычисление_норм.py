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
  - "H1" - $H^1$
  - "H10" - $H^1_0$

- mesh - сетка для вычисления нормы.

```python

```
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

col = st.columns(2)

with col[0]:
    "#### FEniCS"
    with st.echo():
        L2 = norm(u, "L2")
    L2


with col[1]:
    "#### Пример реализации"
    with st.echo():
        from scipy import integrate

        fun = lambda y, x: u(x, y) ** 2
        L2_ = integrate.dblquad(fun, 0, 1, 0, 1)[0] ** 0.5
    L2_

if np.allclose(L2, L2_):
    st.success("Значения достаточно близки")
else:
    st.error("Значения различаются")
