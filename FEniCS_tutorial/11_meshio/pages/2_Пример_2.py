import streamlit as st


r"""
# Пример 1

```python
from fenics import *
import meshio


msh = meshio.read("./Square.msh")
Wri_path = './Dolfin_mesh_functions/'
meshio.write(Wri_path+"mesh.xdmf",
             meshio.Mesh(points=msh.points,
                         cells={"triangle": msh.cells["triangle"]}))

meshio.write(Wri_path+"mf.xdmf",
             meshio.Mesh(points=msh.points,
                         cells={"line": msh.cells["line"]},
                         cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))

mesh = Mesh()
with XDMFFile(Wri_path+"mesh.xdmf") as infile:
    infile.read(mesh)
File(Wri_path+"Dolfin_circle_mesh.pvd").write(mesh)

mvc = MeshValueCollection("size_t", mesh, 1)
with XDMFFile(Wri_path+"mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
File(Wri_path+"Dolfin_circle_facets.pvd").write(mf)

V = FunctionSpace(mesh, 'P', 1)

# Определение граничных условий на основе меток сетки GMSH [Physical Curves: 1, 2, 3, 4]
bc1 = DirichletBC(V, Constant(0.0), mf, 1)
bc2 = DirichletBC(V, Constant(10.0), mf, 2)
bc3 = DirichletBC(V, Constant(6.0), mf, 3)
bc4 = DirichletBC(V, Constant(3.0), mf, 4)

bc = [bc1, bc2, bc3, bc4]

# Определение вариационной задачи
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-15.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Вычисление решения
u = Function(V)
solve(a == L, u, bc)
File("Solution.pvd").write(u)
```

Данный код реализует решение двумерной задачи Дирихле на прямоугольной области. Конечно-элементная сетка создается
внешней библиотекой meshio, которая позволяет читать и записывать файлы сеток в различных форматах. Затем, используя
meshio, создаются два файла формата XDMF, один для сетки и еще один для разметки граничных условий. В последнем файле
разметка граничных условий присваивается ребрам на основе значений меток входящих в cell_data.

Затем, используя встроенный в FEniCS метод XDMFFile, читаются файлы mesh.xdmf и mf.xdmf, после чего создаются объекты
сетки и MeshFunctionSizet (разметка сетки на граничные условия) соответственно.

Далее определяется функциональное пространство и граничные условия. В данном случае, задаются четыре граничных условия
для каждой из сторон прямоугольника. Затем, определяется вариационная задача и ее решение с помощью метода solve.
Решение записывается в файл формата PVD.

Кроме того, этот код записывает различные элементы данных в файлы формата PVD для визуализации с помощью программы
ParaView. В частности, записываются сетка, разметка сетки на граничные условия и решение уравнения.
"""
