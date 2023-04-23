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

gmsh.write("../my_model.msh")

gmsh.finalize()
