import gmsh
from os.path import dirname, join

base_dir = dirname(dirname(__file__))
mesh_dir = join(base_dir, 'mesh')


gmsh.initialize()
gmsh.model.add("my_model")

lc = 0.1  # характеристический размер сетки
p1 = gmsh.model.geo.addPoint(-5, -1, 0, lc)
p2 = gmsh.model.geo.addPoint(5, -1, 0, lc)
p3 = gmsh.model.geo.addPoint(5, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(-5, 1, 0, lc)
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
s = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
gmsh.model.geo.addPlaneSurface([s])

# Создание соответствующей структуры данных Gmsh.
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

for i in range(1):
    gmsh.model.mesh.refine()

# Запуск GUI gmsh
gmsh.fltk.run()

gmsh.option.setNumber("Mesh.MshFileVersion", 2)

gmsh.write(join(mesh_dir, 'my_model.msh'))

#gmsh.write(join(mesh_dir, 'my_model.png'))

gmsh.finalize()


# import gmsh
# from os.path import dirname, join

# base_dir = dirname(dirname(__file__))
# mesh_dir = join(base_dir, 'mesh')
# mesh = Mesh(join(mesh_dir, 'res.xml'))
# #mesh = UnitSquareMesh(200, 40)
# import matplotlib.pyplot as plt

# plot(mesh)
# plt.show()

# print(f'cells {mesh.cells().shape}')
# print(f'vertices {mesh.coordinates().shape}')
# exit(1)