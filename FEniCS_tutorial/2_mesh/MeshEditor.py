from fenics import *
import numpy as np
import matplotlib.pyplot as plt


editor = MeshEditor()

mesh = Mesh()

editor.open(mesh, "triangle", 2, 2)

editor.init_vertices(4)
editor.init_cells(2)

editor.add_vertex(0, Point(0.1, 0.1))
editor.add_vertex(1, np.array([1.0, 0.0]))
editor.add_vertex(2, np.array([0.0, 1.0]))
editor.add_vertex(3, np.array([1.0, 1.0]))

editor.add_cell(0, np.array([0, 1, 3], dtype=np.uintp))
editor.add_cell(1, np.array([0, 2, 3], dtype=np.uintp))

editor.close()

plot(mesh)
plt.show()
