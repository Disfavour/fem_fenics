import streamlit as st


r"""
# Импорт границ и подобластей

Создадим и запишем в файл сетку:

```python
import gmsh

gmsh.initialize()
gmsh.model.add("my_model")

size = 1

lc = 0.1  # характеристический размер сетки

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(size, 0, 0, lc, 2)
gmsh.model.geo.addPoint(size * 2, size, 0, lc, 3)
gmsh.model.geo.addPoint(0, size, 0, lc, 4)

gmsh.model.geo.addLine(1, 2, 5)
gmsh.model.geo.addLine(2, 3, 6)
gmsh.model.geo.addLine(3, 4, 7)
gmsh.model.geo.addLine(4, 1, 8)

gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 20)
gmsh.model.geo.addPlaneSurface([20], 21)

gmsh.model.geo.add_physical_group(1, [5], 31)
gmsh.model.geo.add_physical_group(1, [6], 32)
gmsh.model.geo.add_physical_group(1, [7], 33)
gmsh.model.geo.add_physical_group(1, [8], 34)

gmsh.model.geo.add_physical_group(2, [21], 35)

# Создание соответствующей структуры данных Gmsh.
gmsh.model.geo.synchronize()

# Запуск GUI gmsh
#gmsh.fltk.run()

gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2)

gmsh.write("mshes/square_with_boundaries.msh")

gmsh.finalize()
```

Используем `dolfin-convert`:

```
dolfin-convert square_with_boundaries.msh res.xml
```

Получаем 3 файла:
- res.xml - сетка
- res_physical_region.xml - маркированные подобласти
- res_facet_region.xml - маркированные грани

Пример использования:

```python
from fenics import *
import matplotlib.pyplot as plt

mesh = Mesh("mshes/res.xml")

physical_region = MeshFunction("size_t", mesh, "mshes/res_physical_region.xml")
facet_region = MeshFunction("size_t", mesh, "mshes/res_facet_region.xml")

dx = Measure("dx", domain=mesh, subdomain_data=physical_region)
ds = Measure("ds", domain=mesh, subdomain_data=facet_region)

f = Constant(1)
print(assemble(f * dx(35)))
print(assemble(f * ds(33)))
```
"""
