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

#gmsh.fltk.run()

gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2)
#gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 1)
#gmsh.option.setNumber("Mesh.SaveAll", 0)
#gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)


gmsh.write("mshes/square_with_boundaries.msh")

gmsh.finalize()
