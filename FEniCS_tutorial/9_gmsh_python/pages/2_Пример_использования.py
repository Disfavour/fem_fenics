import streamlit as st


r"""
# Пример использования

Импорт модуля gmsh:

```python
import gmsh
```

Инициализация gmsh и создание нового модельного объекта:

```python
gmsh.initialize()
gmsh.model.add("my_model")
```

Теперь можно создать геометрию вашей сетки, используя команды gmsh. Например, давайте создадим квадрат с вершинами в
точках (0,0), (1,0), (1,1) и (0,1):

```python
lc = 0.1  # характеристический размер сетки
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
s = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
gmsh.model.geo.addPlaneSurface([s])
```

Здесь мы используем функцию `addPoint` для создания четырех вершин квадрата, функцию `addLine` для создания четырех
сторон квадрата, и функцию `addCurveLoop` для объединения сторон в замкнутый контур. Затем мы используем функцию
`addPlaneSurface` для создания плоской поверхности, ограниченной этим контуром.

Для продолжения необходимо создать соответствующую структуру данных gmsh:

```python
gmsh.model.geo.synchronize()
```

Теперь можно задать параметры сетки, например, используя алгоритм Delaunay для автоматической генерации сетки:

```python
gmsh.model.mesh.generate(2)
```

Наконец, можно сохранить сетку в файл, например, в формате msh:

```python
gmsh.write("my_model.msh")
```

В конце необходимо закрыть сессию gmsh:

```python
gmsh.finalize()
```
"""

with st.expander("Программная реализация"):
    r"""
    ```python
import gmsh

gmsh.initialize()
gmsh.model.add("my_model")

lc = 0.1  # характеристический размер сетки
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
s = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
gmsh.model.geo.addPlaneSurface([s])

# Создание соответствующей структуры данных Gmsh.
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("my_model.msh")

gmsh.finalize()
```
    """
