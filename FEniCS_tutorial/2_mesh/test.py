from fenics import *
import matplotlib.pyplot as plt

#File("my_mesh.xml") << UnitSquareMesh(2, 2)

mesh = Mesh("out.xdmf")

plot(mesh)
plt.show()
