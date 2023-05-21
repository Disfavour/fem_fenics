import streamlit as st


r"""
# Конвертация сетки с помощью meshio

```python
import meshio
import numpy as np
import dolfin as df


def gmsh2dolfin(msh_file, mesh_file, line_file=None):

    msh = meshio.gmsh.read(msh_file)

    line_cells = []
    for cell in msh.cells:
        if cell.type == "triangle":
            triangle_cells = cell.data
        elif cell.type == "line":
            if len(line_cells) == 0:
                line_cells = cell.data
            else:
                line_cells = np.vstack([line_cells, cell.data])

    line_data = []
    for key in msh.cell_data_dict["gmsh:physical"].keys():
        if key == "line":
            if len(line_data) == 0:
                line_data = msh.cell_data_dict["gmsh:physical"][key]
            else:
                line_data = np.vstack(
                    [line_data, msh.cell_data_dict["gmsh:physical"][key]]
                )
        elif key == "triangle":
            triangle_data = msh.cell_data_dict["gmsh:physical"][key]

    triangle_mesh = meshio.Mesh(
        points=msh.points,
        cells={"triangle": triangle_cells},
        cell_data={"name_to_read": [triangle_data]},
    )
    line_mesh = meshio.Mesh(
        points=msh.points,
        cells=[("line", line_cells)],
        cell_data={"name_to_read": [line_data]},
    )
    meshio.write(mesh_file, triangle_mesh)
    meshio.xdmf.write(line_file, line_mesh)


def load_mesh(mesh_file):

    mesh = df.Mesh()
    with df.XDMFFile(mesh_file) as infile:
        infile.read(mesh)
        # These are the markers
        ffun = df.MeshFunction("size_t", mesh, 2)
        infile.read(ffun, "name_to_read")

    return mesh, ffun


if __name__ == "__main__":
    msh_file = "biv.msh"
    mesh_file = "mesh.xdmf"
    line_file = "mf.xdmf"

    gmsh2dolfin(msh_file, mesh_file, line_file)
    mesh, ffun = load_mesh(mesh_file)
```

Функция gmsh2dolfin конвертирует файл сетки, созданный в Gmsh (msh_file), в два файла сетки в формате XDMF, которые
могут быть прочитаны библиотекой FEniCS. Один файл содержит треугольные элементы (mesh_file), а другой содержит
линейные элементы (line_file). Эти файлы будут иметь маркеры ячеек, которые могут быть прочитаны FEniCS.

Функция load_mesh загружает сетку из файла XDMF (mesh_file) и возвращает объект сетки (mesh) и метку ячеек (ffun).
Метка ячеек в данном случае относится к маркерам, которые были сохранены в line_file.

Код в основном теле программы использует обе функции для конвертации сетки и загрузки ее в FEniCS.
"""
