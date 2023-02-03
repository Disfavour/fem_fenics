import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Конструктор сеток

## class MeshEditor

Простой редактор для создания сеток.

### Методы класса

- #### open(mesh, type, tdim, gdim)

  Открытие сетки с заданными параметрами.

  ##### Аргументы
  - mesh (Mesh) - сетка.
  - type (str) - тип ячеек сетки.
    - "interval" - для одномерной сетки.
    - "triangle" - для двумерной сетки.
    - "tetrahedron" - для трехмерной сетки.
  - tdim (int) - топологическая размерность ячеек сетки.
  - gdim (int) - геометрическая размерность ячеек сетки.
  
- #### init_vertices(n)

  Указание количества создаваемых вершин.

  ##### Аргументы
  - n (int) - количество вершин.
  
- #### init_cells(n)

  Указание количества создаваемых ячеек.

  ##### Аргументы
  - n (int) - количество ячеек.

- #### add_vertex(index, x)
  
  Добавление вершины по заданным координатам.
  
  ##### Аргументы
  - index (int) - номер создаваемой вершины.
  - x (numpy.array(float) или Point) - координаты точки.

- #### add_cell(c, v)
  
  Добавление ячейки с заданными вершинами
  
  ##### Аргументы
  - c (int) - номер создаваемой ячейки.
  - v (numpy.array(int)) - номера вершин, которые образуют ячейку.

- #### close()

  Закрытие сетки.
"""

with st.expander("Пример использования MeshEditor"):
    with st.echo():
        from fenics import *
        import numpy as np
        import matplotlib.pyplot as plt

        # Создание объекта класса MeshEditor
        editor = MeshEditor()

        # Создание пустой сетки
        mesh = Mesh()

        editor.open(mesh, "triangle", 2, 2)

        editor.init_vertices(4)
        editor.init_cells(2)

        editor.add_vertex(0, Point(0.0, 0.0))
        editor.add_vertex(1, np.array([1.0, 0.0]))
        editor.add_vertex(2, np.array([0.0, 1.0]))
        editor.add_vertex(3, np.array([1.0, 1.0]))

        editor.add_cell(0, np.array([0, 1, 3], dtype=np.uintp))
        editor.add_cell(1, np.array([0, 2, 3], dtype=np.uintp))

        editor.close()

        plot(mesh)

    col = st.columns(2)

    with col[0]:
        st.pyplot(plt.gcf())
        plt.clf()

with st.expander("Удаление ячейки с помощью MeshEditor"):
    with st.echo():
        from fenics import *
        import numpy as np
        import matplotlib.pyplot as plt

        mesh = UnitIntervalMesh(5)

        old_cells = mesh.cells()
        old_coors = mesh.coordinates()

        new_coors = np.delete(old_coors, (3), axis=0)
        new_cells = np.delete(old_cells, (-1), axis=0)

        new_mesh = Mesh()
        editor = MeshEditor()
        editor.open(new_mesh, "interval", 1, 1)

        editor.init_vertices(len(new_coors))
        editor.init_cells(len(new_cells))

        vert_id = 0
        for vert in new_coors:
            editor.add_vertex(vert_id, vert)
            vert_id += 1

        cell_id = 0
        for c in range(len(new_cells)):
            editor.add_cell(cell_id, new_cells[c])
            cell_id += 1

        editor.close()

        plot(new_mesh)

    st.pyplot(plt.gcf())
    plt.clf()
