import streamlit as st
from PIL import Image
import os.path


domain112 = Image.open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "images", "domain112.png"))


r"""
# Множественные условия Дирихле

## Краевая задача

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- & \div \grad u = f \quad & \bold{x} & \in \Omega,
\\[0.5 cm]
& u = u_L \quad & \bold{x} & \in \Gamma_D^L,
\\[0.5 cm]
& u = u_R \quad & \bold{x} & \in \Gamma_D^R,
\\[0.5 cm]
& \frac {\partial u} {\partial \bm{n}} = g \quad & \bold{x} & \in \Gamma_N.
\end{aligned}
$$


## Вариационная формулировка

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\def \phi {\varphi}
\def \f {\bold{f}}
\def \g {\bold{g}}

\begin{aligned}

- &\div \grad u = f
\\[0.5 cm]
- &\int \limits_\Omega v \div \grad u \ dx = \int \limits_\Omega f v \ dx
\\[0.5 cm]
- &\int \limits_\Omega \div (v \grad u) \ dx + \int \limits_\Omega \grad u \cdot \grad v \ dx
= \int \limits_\Omega f v \ dx
\\[0.5 cm]
&\int \limits_\Omega \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_\Omega \div (v \grad u) \ dx
\\[1.0 cm]
&\text{По теореме о диверенции:}
\\[0.5 cm]
&\int \limits_\Omega \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_{\partial \Omega} \frac {\partial u} {\partial \bm{n}} v \ ds
\\[1.0 cm]
&\text{Так как } v = 0 \text{ на } \Gamma_D:
\\[0.5 cm]
&\int \limits_\Omega \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_{\Gamma_N} \frac {\partial u} {\partial \bm{n}} v \ ds
\\[0.5 cm]
&\int \limits_\Omega \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_{\Gamma_N} g v \ ds
\\[1.0 cm]
&\text{В стандартной записи:}
\\[0.5 cm]
&a(u, v) = \int \limits_\Omega \grad u \cdot \grad v \ dx
\\[0.5 cm]
&L(v) = \int \limits_\Omega f v \ dx + \int \limits_{\Gamma_N} g v \ ds

\end{aligned}
$$


## Тестовая двумерная задача
"""

with st.columns(3)[1]:
    st.image(domain112)

r"""
$$
\begin{aligned}
u_e &= 1 + x_1^2 + 2 x_2^2,
\\[0.5 cm]
f(x_1, x_2) &= -6,
\\[0.5 cm]
g(x_1, x_2) &= 
\begin{cases}
0, &x_2 = 0 \\
4, &x_2 = 1
\end{cases},
\\[0.5 cm]
u_L(x_1, x_2) &= 1 + 2 x_2^2,
\\[0.5 cm]
u_R(x_1, x_2) &= 2 + 2 x_2^2.
\end{aligned}
$$

Можно определить $g$ как функцию на области так, чтобы значения на границах были корректными:

$$
g(x_1, x_2) = 4 x_2.
$$


## Необходимые изменения относительно базовой версии

Для $\Gamma_D^L$ определяем граничное условие:

```python
tol = 1E-14

u_L = Expression("1 + 2*x[1]*x[1]", degree=2)

def boundary_L(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol)
    
bc_L = DirichletBC(V, u_L, boundary_L)
```

Для $\Gamma_D^R$ аналогично:

```python
u_L = Expression("2 + 2*x[1]*x[1]", degree=2)

def boundary_L(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 1, tol)
    
bc_L = DirichletBC(V, u_L, boundary_L)
```

Собираем граничные условия в список, который можно передать функции solve для вычисления решения:

```python
bcs = [bc_L, bc_R]
...
solve(a == L, u, bcs)
```
"""

with st.expander("Программа"):
    """
    ```python
    from fenics import *
    import matplotlib.pyplot as plt
    
    
    mesh = UnitSquareMesh(8, 8)
    
    V = FunctionSpace(mesh, "P", 1)
    
    tol = 1E-14
    
    u_L = Expression("1 + 2*x[1]*x[1]", degree=2)
    
    
    def boundary_L(x, on_boundary):
        return on_boundary and near(x[0], 0, tol)
    
    
    bc_L = DirichletBC(V, u_L, boundary_L)
    
    u_R = Expression("2 + 2*x[1]*x[1]", degree=2)
    
    
    def boundary_R(x, on_boundary):
        return on_boundary and near(x[0], 1, tol)
    
    
    bc_R = DirichletBC(V, u_R, boundary_R)
    
    bcs = [bc_L, bc_R]
    
    # exact solution
    u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v))*dx
    g = Expression("-4*x[1]", degree=1)
    L = f*v*dx - g*v*ds
    
    u = Function(V)
    solve(a == L, u, bcs)
    
    plt.colorbar(plot(u))
    plot(mesh)
    
    plt.show()
    ```
    """


r"""
## Параметрические расчеты

"""

col = st.columns([2, 3])

with col[0]:
    nx1 = st.slider("nx1", 1, 100, 10)
    nx2 = st.slider("nx2", 1, 100, 10)
    diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"), 2)
    degree = st.slider("degree", 1, 3, 1)
    show_mesh = st.checkbox('визуализация сетки')


from fenics import *
import numpy as np
import matplotlib.pyplot as plt


mesh = UnitSquareMesh(nx1, nx2, diagonal)

V = FunctionSpace(mesh, "P", degree)

tol = 1E-14

u_L = Expression("1 + 2*x[1]*x[1]", degree=2)


def boundary_L(x, on_boundary):
    return on_boundary and near(x[0], 0, tol)


bc_L = DirichletBC(V, u_L, boundary_L)

u_R = Expression("2 + 2*x[1]*x[1]", degree=2)


def boundary_R(x, on_boundary):
    return on_boundary and near(x[0], 1, tol)


bc_R = DirichletBC(V, u_R, boundary_R)

bcs = [bc_L, bc_R]

# exact solution
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
g = Expression("-4*x[1]", degree=1)
L = f*v*dx - g*v*ds

u = Function(V)
solve(a == L, u, bcs)

plt.colorbar(plot(u))
if show_mesh:
    plot(mesh)

error_L2 = errornorm(u_D, u, "L2")

with col[1]:
    st.pyplot(plt.gcf())

with st.columns([1, 2, 1])[1]:
    f"""
    Норма ошибки $L^2 = {error_L2}$
    """

