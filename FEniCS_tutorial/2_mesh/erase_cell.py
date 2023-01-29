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
plt.show()