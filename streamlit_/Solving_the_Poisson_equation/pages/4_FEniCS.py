import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# FEniCS

### Подключение библиотеки

```python
from fenics import *
```
"""

r"""
### Создание сетки

#### Point(x1=0.0, x2=0.0, x3=0.0)

Создание точки $(x_1, x_2, x_3)$.
"""

r"""
#### IntervalMesh(n, a, b)

Создание равномерной сетки на отрезке $[a, b]$.

##### Параметры

- n (int) - количество ячеек.

- a (float) - левая граница отрезка.

- b (float) - правая граница отрезка.
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        n = st.slider("n", 1, 100, 10)

    with col2:
        a, b = st.slider("a, b", -10.0, 10.0, (0.0, 1.0))

    plot(IntervalMesh(n, a, b))
    st.pyplot(plt.gcf())
    plt.clf()

r"""
#### UnitIntervalMesh(n)

Аналог IntervalMesh, но только для отрезка $[0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitIntervalMesh(n)
```

```python
IntervalMesh(n, 0, 1)
```
"""

r"""
#### RectangleMesh(p1, p2, nx1, nx2, diagonal="right")

Создание треугольной равномерной сетки на прямоугольнике.

##### Параметры

- p1 (Point) - точка с минимальными значениями $x_1, x_2$ прямоугольника.

- p2 (Point) - точка с максимальными значениями $x_1, x_2$ прямоугольника.

- nx1 (int) - количество ячеек по $x_1$.

- nx2 (int) - количество ячеек по $x_2$.

- diagonal (str) - строка, указывающая направление диагоналей ячеек.
  
  Может принимать следующие значения:

    - "left"
    - "right"
    - "right/left"
    - "left/right"
    - "crossed"
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        "##### p1"
        p1 = Point(st.slider("x1", -10.0, 10.0, 0.0, key=11), st.slider("x2", -10.0, 10.0, 0.0, key=12))

    with col2:
        "##### p2"
        p2 = Point(st.slider("x1", -10.0, 10.0, 1.0, key=13), st.slider("x2", -10.0, 10.0, 1.0, key=14))

    col1, col2, col3 = st.columns(3)

    with col1:
        nx1 = st.slider("nx1", 1, 100, 10, key=15)

    with col2:
        nx2 = st.slider("nx2", 1, 100, 10, key=16)

    with col3:
        diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"))

    plot(RectangleMesh(p1, p2, nx1, nx2, diagonal))
    st.pyplot(plt.gcf())
    plt.clf()

r"""
#### UnitSquareMesh(x1n, x2n, diagonal="right")

Аналог RectangleMesh, но только для квадрата $[0, 1] \times [0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitSquareMesh(x1n, x2n, diagonal)
```

```python
RectangleMesh(Point(0, 0), Point(1, 1), x1n, x2n, diagonal)
```
"""

r"""
#### BoxMesh(p1, p2, nx1, nx2, nx3)

Создание равномерной сетки из тетраэдров на прямоугольном параллелепипеде.

##### Параметры

- p1 (Point) - точка с минимальными значениями $x_1, x_2, x_3$ прямоугольного параллелепипеда.
  
- p2 (Point) - точка с максимальными значениями $x_1, x_2, x_3$ прямоугольного параллелепипеда.

- nx1 (int) - количество ячеек по $x_1$.

- nx2 (int) - количество ячеек по $x_2$.

- nx3 (int) - количество ячеек по $x_3$.
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        "##### p1"
        p1 = Point(
            st.slider("x1", -10.0, 10.0, 0.0, key=21),
            st.slider("x2", -10.0, 10.0, 0.0, key=22),
            st.slider("x3", -10.0, 10.0, 0.0, key=23),
        )

    with col2:
        "##### p2"
        p2 = Point(
            st.slider("x1", -10.0, 10.0, 1.0, key=24),
            st.slider("x2", -10.0, 10.0, 1.0, key=25),
            st.slider("x3", -10.0, 10.0, 1.0, key=26),
        )

    col1, col2, col3 = st.columns(3)

    with col1:
        nx1 = st.slider("nx1", 1, 100, 10, key=27)

    with col2:
        nx2 = st.slider("nx2", 1, 100, 10, key=28)

    with col3:
        nx3 = st.slider("nx3", 1, 100, 10, key=29)

    plot(BoxMesh(p1, p2, nx1, nx2, nx3))
    st.pyplot(plt.gcf())
    plt.clf()

r"""
#### UnitCubeMesh(nx1, nx2, nx3)

Аналог BoxMesh, но только для куба $[0, 1] \times [0, 1] \times [0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitCubeMesh(nx1, nx2, nx3)
```

```python
BoxMesh(Point(0, 0, 0), Point(1, 1, 1), nx1, nx2, nx3)
```
"""

r"""
#### Сетка из файла

Создание сетки из файла в формате DOLPHIN XML:

```python
Mesh("file.xml")
```
"""

r"""
### Задание конечно-элементного функционального пространства

Для создания функционального пространства используется следующая функция:

```python
V = FunctionSpace(mesh, "P", 1)
```

Второй аргумент "P" определяет тип элемента. Тип элемента P подразумевает стандартное семейство элементов Лагранжа.
Третий аргумент определяет степень конечного элемента.

### Задание пробных и тестовых функций

В FEniCS не указываются граничные условия как часть функционального пространства, поэтому достаточно работать с одним 
общим пространством $V$ как для пробных, так и для тестовых функций в программе.

```python
u = TrialFunction(V)
v = TestFunction(V)
```

### Задание граничных условий

```python
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, g, boundary)
```

Функция boundary определяет, какие точки принадлежат той части границы, к которой должно быть применено граничное 
условие.

Аргумент on_boundary имеет значение True, если $\bm x$ находится на физической границе сетки, поэтому можно просто
вернуть его.

Функция boundary будет вызываться для каждой дискретной точки в сетке, что означает, что при желании мы можем 
определить границы, где $u$ также известно внутри области.

При определении Expression второй аргумент степень является параметром, который определяет, как выражение должно
обрабатываться в вычислениях. Для каждого локального элемента FEniCS интерполирует выражение в пространство конечных
элементов заданной степени. Чтобы получить оптимальную точность вычислений, обычно хорошим выбором является
использование той же степени, что и для пространства $V$, которое используется для пробных и тестовых функций. Однако,
если выражение используется для представления точного решения, которое используется для оценки точности вычисленного
решения, для выражения должна использоваться более высокая степень (на одну или две степени выше).

### Задание правой части

```python
f = Expression(’-6’, degree=0)
```

или

```python
f = Constant(-6)
```

### Задание вариационной задачи

```python
a = inner(grad(u), grad(v))*dx
L = f*v*dx
```

### Формирование и решение СЛАУ

```python
u = Function(V)
solve(a == L, u, bc)
```

### Визуализация

```python
plot(u)
plot(mesh)
plt.show()
```

### Экспорт решения в формате VTK

```python
vtkfile = File("poisson/solution.pvd")
vtkfile << u
```

### Вычисление ошибки

Для проверки точности решения можно вычислить $L^2$ норму ошибки.

```python
error_L2 = errornorm(u_D, u, "L2")
```

Вычисления максимального значения ошибки в вершинах:

```python
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
```
"""

"""
### Программа

```python
from fenics import *
import matplotlib.pyplot as plt
import numpy as np


# Создание сетки и функционального пространства
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)

# Граничные условия
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, g, boundary)

# Вариационная постановка
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
l = f*v*dx

# Вычисление
u = Function(V)
solve(a == l, u, bc)

# Визуализация сетки и решения
plt.colorbar(plot(u))
plot(mesh)
plt.show()

# Экспорт решения в формате VTK
vtkfile = File("poisson/solution.pvd")
vtkfile << u

# Вычисление ошибки
error_L2 = errornorm(u_D, u, "L2")

vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
```
"""
