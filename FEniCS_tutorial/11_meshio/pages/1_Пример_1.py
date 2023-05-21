import streamlit as st


r"""
# Пример 1

Вот простой пример использования meshio для чтения файла сетки и использования его с FEniCS:

```python
import meshio
from dolfin import *
import matplotlib.pyplot as plt

# Читаем файл сетки используя meshio
mesh = meshio.read("mesh.xml")

# Конвертируем данные сетки в формат FEniCS
vertices = mesh.points
cells = {"triangle": mesh.cells["triangle"]}

# Создаем объект сетки FEniCS
mesh_fenics = Mesh(MeshValueCollection("size_t", 2, mesh.cells["triangle"]), vertices)

# Определяем функциональное пространство
V = FunctionSpace(mesh_fenics, "P", 1)

# Определяем граничные условия
def boundary(x, on_boundary):
    return on_boundary

u_D = Constant(0.0)
bc = DirichletBC(V, u_D, boundary)

# Определяем вариационную задачу
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# Решаем вариационную задачу
u = Function(V)
solve(a == L, u, bc)

# Визуализация решения
plot(u)
plt.show()
```

Здесь мы сначала используем meshio для чтения файла сетки в формате XML, затем конвертируем данные сетки в формат
FEniCS, создаем объект сетки FEniCS и определяем функциональное пространство. Затем мы определяем граничные условия,
вариационную задачу и решаем ее, используя функцию solve. Наконец, мы визуализируем решение с помощью функции plot.
"""
