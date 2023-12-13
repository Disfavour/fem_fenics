from fenics import *
from os.path import dirname, join

base_dir = dirname(dirname(__file__))
mesh_dir = join(base_dir, 'mesh')
mesh = Mesh(join(mesh_dir, 'res.xml'))
#mesh = UnitSquareMesh(200, 40)
import matplotlib.pyplot as plt

#plt.figure(figsize=(6.4, 3.6), dpi=600, tight_layout=True)
plot(mesh)
plt.show()

# plt.savefig(join(mesh_dir, f'mesh.png'), dpi=600, tight_layout=True)
    
# plt.close()

print(f'cells {mesh.cells().shape}')
print(f'vertices {mesh.coordinates().shape}')
