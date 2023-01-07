import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Решение с помощью FEniCS

### Подключение библиотеки

```python
from fenics import *
```

### Создание сетки

Для создания сетки n $\times$ m используется следующая функция:

```python
mesh = UnitSquareMesh(n, m)
```

Каждая прямоугольная ячейка состоит из 2-х треугольников.
"""


if st.checkbox("Визцализация сетки"):
    n = st.slider("Размерность сетки n", 1, 100, 10)
    m = st.slider("Размерность сетки m", 1, 100, 10)
    plot(UnitSquareMesh(n, m))
    st.pyplot(plt.gcf())


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
