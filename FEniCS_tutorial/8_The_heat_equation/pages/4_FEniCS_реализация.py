import streamlit as st


r"""
# FEniCS реализация

Главные граничные условия реализуются обычным образом:

```python
bc = DirichletBC(V, u_D, "on_boundary")
```

Будем использовать `y` для $y^{n+1}$ на новом временном шаге и `y_n` для $y^n$.

Начальное значение `y_n` может быть вычислено с помощью проекции `project`, что является аналогом решения
соответствующей вариационной задачи:

```python
y_n = project(Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))", degree=2, a=5), V)
```

Мы можем либо определить $a$ или $L$, либо мы можем просто определить $F$ и попросить FEniCS выяснить, какие члены
должны входить в билинейную форму $a$, а какие - в линейную форму $L$.
Последнее удобно, особенно в более сложных задачах.

```python
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)
a = tau*dot(grad(y), grad(v))*dx + y*v*dx
L = y_n*v*dx + tau*f*v*dx
```

Шаги во временном в цикле:

```python
y = Function(V)
t = 0
for n in range(num_steps):
    # Обновляем время
    t += tau
    # Вызываем решатель
    solve(a == L, y, bc)
    # Обновляем предыдущее решение
    y_n.assign(y)
```

На последнем шаге цикла с шагом по времени мы присваиваем значения переменной `y` (новое вычисленное решение)
переменной `y_n`, содержащей значения на предыдущем шаге по времени. Это должно быть сделано с помощью метода
`assign`. Если мы вместо этого попытаемся сделать `y_n = y`, мы установим переменную `y_n` такой же переменной, как `y`,
что не является тем, что мы хотим.
(Нам нужны две переменные, одна для значений на предыдущем временном шаге и одна для значений на текущем
временном шаге.)
"""

with st.expander("Программная реализация"):
    r"""
```python
from fenics import *
import matplotlib.pyplot as plt


T = 2.0
num_steps = 50
tau = T / num_steps

nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)

V = FunctionSpace(mesh, "P", 1)

bc = DirichletBC(V, Constant(0), "on_boundary")

y_n = project(Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))", degree=2, a=5), V)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = tau*dot(grad(u), grad(v))*dx + u*v*dx
L = y_n*v*dx + tau*f*v*dx

y = Function(V)
t = 0
for n in range(num_steps):
    # Обновляем время
    t += tau
    # Считаем u на следующем шаге
    solve(a == L, y, bc)

    plot(y)
    plt.show()

    # Обновляем результат предыдущего шага
    y_n.assign(y)
```
    """
