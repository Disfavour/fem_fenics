import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# FEniCS

### Подключение библиотеки

```python
from fenics import *
```

### Создание сетки

Для создания сетки $n \times m$ ячеек на области $[0, 1] \times [0, 1]$ используется следующая функция:

```python
mesh = UnitSquareMesh(n, m)
```

Каждая ячейка состоит из двух треугольников.
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        n = st.slider("n", 1, 100, 10)

    with col2:
        m = st.slider("m", 1, 100, 10)

    plot(UnitSquareMesh(n, m))
    st.pyplot(plt.gcf())
    plt.clf()

r"""
### Задание функционального пространства

Для создания функционального пространства используется следующая функция:

```python
V = FunctionSpace(mesh, t, n)
```

- mesh - сетка.
- t - тип конечных элементов ("P" - лагранжевы конечные элементы).
- n - степень конечных элементов.

### Задание пробных и тестовых функций

В FEniCS не указываются граничные условия как часть функционального пространства, поэтому достаточно работать с одним 
общим пространством $V$ как для пробных, так и для тестовых функций в программе.

```python
u = TrialFunction(V)
v = TestFunction(V)
```

### Задание граничных условий

Для задания граничных условий используется следующая функция:

```python
bc = DirichletBC(V, g, boundary)
```

- V - функциональное пространство.
- g - Expression, определяющее значения на границе.
- boundary - функция, определяещая какие точки принадлежат той части границы, к которой должно быть применено граничное 
  условие.

g задается следующим образом:

```python
g = Expression(formula, degree=n)
```

- formula - строка, содержащая математическое выражение в синтаксисе С++. Вместо x, y, z ипользуются x[0], x[1], x[2].
  Например, "1 + x[0] * x[0] + 2 * x[1] * x[1]".
- n - степень интерполированного выражения в пространстве конечных элементов. Оптимальным является использование той же
  степени, что и для пространства V.
  
Функция boundary может быть задана следующим образом:

```python
def boundary(x, on_boundary):
    return on_boundary
```

- x - точка, для которой необходимо определить принадлежность границе.
- on_boundary - вспомогательный аргумент, который принимает значение `True`, если точка x на границе сетки, и `False`
  иначе.

Функция должна возвращать `True`, если точка x лежит на границе, и `False` в противном случае.

Функция boundary будет вызываться для каждой дискретной точки в сетке, поэтому можно  определить границы, где $u$
также известно внутри области.

### Задание правой части

Правая часть задается с помощью Expression:

```python
f = Expression(formula, degree=n)
```

Когда правая часть константа (const), более эффективно использовать Constant:

```python
f = Constant(const)
```

### Задание вариационной задачи

Теперь есть все необходимое для определения вариационной задачи:

```python
a = inner(grad(u), grad(v)) * dx
l = f * v * dx
```

### Формирование и решение СЛАУ

Определив конечно-элементную вирационную задачу и граничные условия, можно получить решение:

```python
u = Function(V)
solve(a == l, u, bc)
```

Function - объект, представляющий решение.

### Визуализация

Для визуализации решения и сетки можно использовать фукнцию plot:

```python
plot(u)
plot(mesh)
```

Для получения интеркативного изображения необходим следующий вызов:

```python
plt.show()
```

### Экспорт решения в формате VTK

Решение можно экспортировать следующим образом:

```python
vtkfile = File("poisson/solution.pvd")
vtkfile << u
```

Например, для использования ParaView для визуализации.

### Вычисление ошибки

Для проверки точности решения можно вычислить $L^2$ норму ошибки.

```python
error_L2 = errornorm(u_D, u, "L2")
```

Вычисления максимального значения ошибки в вершинах:

```python
vertex_values_g = g.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_g - vertex_values_u))
```

Вычисляются значения $g$ и $u$ во всех вершинах, затем результаты вычитаются по модулю.
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

# Вариационная задача
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(grad(u), grad(v))*dx
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
