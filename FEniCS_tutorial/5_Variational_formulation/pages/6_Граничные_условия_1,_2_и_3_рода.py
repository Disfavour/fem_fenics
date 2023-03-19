import streamlit as st
from PIL import Image
import os.path


domain123 = Image.open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "images", "domain123.png"))


r"""
# Граничные условия Дирихле и Неймана

## Краевая задача

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- & \div (\kappa \grad u) = f \quad & \bold{x} & \in \Omega,
\\[0.5 cm]
& u = u_D^i \quad & \bold{x} & \in \Gamma_D^i, \quad i = 0, 1, \dots
\\[0.5 cm]
& \kappa \frac {\partial u} {\partial \bm{n}} = g_i \quad & \bold{x} & \in \Gamma_N^i, \quad i = 0, 1, \dots
\\[0.5 cm]
& \kappa \frac {\partial u} {\partial \bm{n}} + \alpha u = q_i \quad & \bold{x} & \in \Gamma_R^i, \quad i = 0, 1, \dots
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

- & \div (\kappa \grad u) = f
\\[0.5 cm]
- &\int \limits_\Omega v \div (\kappa \grad u) \ dx = \int \limits_\Omega f v \ dx
\\[0.5 cm]
- &\int \limits_\Omega \div (v \kappa \grad u) \ dx + \int \limits_\Omega \kappa \grad u \cdot \grad v \ dx
= \int \limits_\Omega f v \ dx
\\[0.5 cm]
&\int \limits_\Omega \kappa \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_\Omega \div (v \kappa \grad u) \ dx
\\[1.0 cm]
&\text{По теореме о диверенции:}
\\[0.5 cm]
&\int \limits_\Omega \kappa \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \int \limits_{\partial \Omega} \kappa \frac {\partial u} {\partial \bm{n}} v \ ds
\\[1.0 cm]
&\text{Так как } v = 0 \text{ на } \Gamma_D^i:
\\[0.5 cm]
&\int \limits_\Omega \kappa \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \sum_i \int \limits_{\Gamma_N^i} \kappa \frac {\partial u} {\partial \bm{n}} v \ ds
+ \sum_i \int \limits_{\Gamma_R^i} \kappa \frac {\partial u} {\partial \bm{n}} v \ ds
\\[0.5 cm]
&\int \limits_\Omega \kappa \grad u \cdot \grad v \ dx = \int \limits_\Omega f v \ dx
+ \sum_i \int \limits_{\Gamma_N^i} g_i v \ ds
+ \sum_i \int \limits_{\Gamma_R^i} (q_i - \alpha u) v \ ds
\\[1.0 cm]
&F = \int \limits_\Omega \kappa \grad u \cdot \grad v \ dx - \int \limits_\Omega f v \ dx
- \sum_i \int \limits_{\Gamma_N^i} g_i v \ ds
- \sum_i \int \limits_{\Gamma_R^i} (q_i - \alpha u) v \ ds = 0
\\[0.5 cm]
&F = \int \limits_\Omega \kappa \grad u \cdot \grad v \ dx
- \int \limits_\Omega f v \ dx
- \sum_i \int \limits_{\Gamma_N^i} g_i v \ ds
- \sum_i \int \limits_{\Gamma_R^i} q_i v \ ds
+ \sum_i \int \limits_{\Gamma_R^i} \alpha u v \ ds = 0
\\[1.0 cm]
&\text{В стандартной записи:}
\\[0.5 cm]
&a(u, v) = \int \limits_\Omega \kappa \grad u \cdot \grad v \ dx + \sum_i \int \limits_{\Gamma_R^i} \alpha u v \ ds
\\[0.5 cm]
&L(v) = \int \limits_\Omega f v \ dx
+ \sum_i \int \limits_{\Gamma_N^i} g_i v \ ds
+ \sum_i \int \limits_{\Gamma_R^i} q_i v \ ds

\end{aligned}
$$


## Тестовая двумерная задача
"""

with st.columns(3)[1]:
    st.image(domain123)

r"""
$$
\begin{aligned}
u_e &= 1 + x_1^2 + 2 x_2^2,
\\[0.5 cm]
\kappa &= 1,
\\[0.5 cm]
f(x_1, x_2) &= -6,
\\[0.5 cm]
g(x_1, x_2) &= 4 x_2 = 4,
\\[0.5 cm]
u_D(x_1, x_2) &= 1 + x_1^2 + 2 x_2^2,
\\[1.0 cm]
\kappa \frac {\partial u} {\partial \bm{n}} + \alpha u &= q,
\\[0.5 cm]
\alpha u &= q,
\\[0.5 cm]
\alpha (u - \frac {q} {\alpha}) &= 0,
\\[0.5 cm]
\alpha &= 1000,
\\[0.5 cm]
\frac {q} {\alpha} &= u_e.
\end{aligned}
$$


## Необходимые изменения относительно базовой версии

Из выражения $F$ можно получить $a(u, v)$ и $L(v)$ с помощью функций:

```python
a = lhs(F)
L = rhs(F)
```

Промаркеруем соответстующие части границы:

```python
boundary_markers = FacetFunction("size_t", mesh)

class BoundaryX0(SubDomain):
    tol = 1E-14
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0, tol)
        
bx0 = BoundaryX0()
bx0.mark(boundary_markers, 0)
```

Аналогично для каждой границы.

Для использования интегралов по частям границы $ds(i)$ нужно переопределить меру $ds$ в терминах наших граничных
маркеров:

```python
ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
```

Аналогично, для использования интегралов по разным частям области нужно переопределить $dx$:

```python
dx = Measure("dx", domain=mesh, subdomain_data=domain_markers)
```

Предположим, у нас есть условие Робина со значениями r и s на подобласти R и условие Неймана со значением g
на подобласти N. Вариационная форма может быть записана:

```python
a = kappa*dot(grad(u), grad(v))*dx + alpha*u*v*ds(R)
L = f*v*dx + g*v*ds(N) + q*v*ds(R)
```
"""

with st.expander("Программа"):
    """
    ```python
    from fenics import *
    import matplotlib.pyplot as plt
    
    
    tol = 1E-14
    kappa = 1
    alpha = 1000
    
    
    mesh = UnitSquareMesh(8, 8)
    
    V = FunctionSpace(mesh, "P", 1)
    
    boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    
    bx0 = CompiledSubDomain("on_boundary && near(x[0], 0, tol)", tol=tol)
    bx1 = CompiledSubDomain("on_boundary && near(x[0], 1, tol)", tol=tol)
    by0 = CompiledSubDomain("on_boundary && near(x[1], 0, tol)", tol=tol)
    by1 = CompiledSubDomain("on_boundary && near(x[1], 1, tol)", tol=tol)
    
    bx0.mark(boundary_markers, 0)
    bx1.mark(boundary_markers, 1)
    by0.mark(boundary_markers, 2)
    by1.mark(boundary_markers, 3)
    
    u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    
    bcs = [
        DirichletBC(V, u_D, bx0),
        DirichletBC(V, u_D, bx1),
    ]
    
    ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    
    g = Constant(4)
    
    q = Expression(f"(1 + x[0]*x[0] + 2*x[1]*x[1]) * {alpha}", degree=2)
    
    F = kappa*dot(grad(u), grad(v))*dx + alpha * u * v * ds(2) - f*v*dx - g*v*ds(3) - q*v*ds(2)
    a = lhs(F)
    L = rhs(F)
    
    # a = kappa*dot(grad(u), grad(v))*dx + alpha * u * v * ds(2)
    # L = f*v*dx + g*v*ds(3) + q*v*ds(2)
    
    u = Function(V)
    solve(a == L, u, bcs)
    
    plt.colorbar(plot(u))
    plot(mesh)
    
    error_L2 = errornorm(u_D, u, "L2")
    
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
import matplotlib.pyplot as plt

tol = 1E-14
kappa = 1
alpha = 1000

mesh = UnitSquareMesh(nx1, nx2, diagonal)

V = FunctionSpace(mesh, "P", degree)

boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)

bx0 = CompiledSubDomain("on_boundary && near(x[0], 0, tol)", tol=tol)
bx1 = CompiledSubDomain("on_boundary && near(x[0], 1, tol)", tol=tol)
by0 = CompiledSubDomain("on_boundary && near(x[1], 0, tol)", tol=tol)
by1 = CompiledSubDomain("on_boundary && near(x[1], 1, tol)", tol=tol)

bx0.mark(boundary_markers, 0)
bx1.mark(boundary_markers, 1)
by0.mark(boundary_markers, 2)
by1.mark(boundary_markers, 3)

u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

bcs = [
    DirichletBC(V, u_D, bx0),
    DirichletBC(V, u_D, bx1),
]

ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)

g = Constant(4)

q = Expression(f"(1 + x[0]*x[0] + 2*x[1]*x[1]) * {alpha}", degree=2)

F = kappa * dot(grad(u), grad(v)) * dx + alpha * u * v * ds(2) - f * v * dx - g * v * ds(3) - q * v * ds(2)
a = lhs(F)
L = rhs(F)

# a = kappa*dot(grad(u), grad(v))*dx + alpha * u * v * ds(2)
# L = f*v*dx + g*v*ds(3) + q*v*ds(2)

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

