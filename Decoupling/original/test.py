from dolfin import *

# u.compute_vertex_values()
# mesh.coordinates()
# u.vector().get_local()
# mesh.geometric_dimension()
# u.vector()[:]

n = 3
mesh = UnitSquareMesh(n, n)

n = FacetNormal(mesh)

g = Constant((100, 100))  # Expression(('1', '1'), degree=1)


# P1
VP = VectorFunctionSpace(mesh, 'P', 1)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

subdomain_left = CompiledSubDomain('on_boundary && near(x[0], 0)')
subdomain_right = CompiledSubDomain('on_boundary && near(x[0], 1)')
subdomain_top = CompiledSubDomain('on_boundary && near(x[1], 1)')
subdomain_bot = CompiledSubDomain('on_boundary && near(x[1], 0)')

subdomain_left.mark(boundaries, 1)
subdomain_right.mark(boundaries, 2)
subdomain_top.mark(boundaries, 3)
subdomain_bot.mark(boundaries, 4)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

bcs = [DirichletBC(VP.sub(0), Constant(0), boundaries, 1),
       DirichletBC(VP.sub(0), Constant(0), boundaries, 2),
       DirichletBC(VP.sub(1), Constant(0), boundaries, 3),
       DirichletBC(VP.sub(1), Constant(0), boundaries, 4)]

# bcs = [DirichletBC(VP.sub(0), Constant(0), 'on_boundary && near(x[0], 0)'),
#        DirichletBC(VP.sub(0), Constant(0), 'on_boundary && near(x[0], 1)'),
#        DirichletBC(VP.sub(1), Constant(0), 'on_boundary && near(x[1], 0)'),
#        DirichletBC(VP.sub(1), Constant(0), 'on_boundary && near(x[1], 1)')]

u = project(g, VP)

print(assemble(dot(u, n) * ds))

bcs[0].apply(u.vector())
bcs[1].apply(u.vector())

print(assemble(dot(u, n) * ds))

#print(u.compute_vertex_values())

# RT
RT = FunctionSpace(mesh, 'RT', 1)

# boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
# subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], 0) || near(x[0], 1))")
# subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], 0) || near(x[1], 1))")
# subdomain_x.mark(boundaries, 1)
# subdomain_y.mark(boundaries, 2)
#

# bcs = [DirichletBC(VP, Constant((0, 0)), 'on_boundary && near(x[0], 0)'),
#        DirichletBC(VP, Constant((0, 0)), 'on_boundary && near(x[0], 1)'),
#        DirichletBC(VP, Constant((0, 0)), 'on_boundary && near(x[1], 0)'),
#        DirichletBC(VP, Constant((0, 0)), 'on_boundary && near(x[1], 1)')]

bcs = [DirichletBC(RT, Constant((0, 0)), boundaries, 1),
       DirichletBC(RT, Constant((0, 0)), boundaries, 2),
       DirichletBC(RT, Constant((0, 0)), boundaries, 3),
       DirichletBC(RT, Constant((0, 0)), boundaries, 4)]

#bc = DirichletBC(RT, Constant((0, 0)), 'on_boundary')

u = project(g, RT)
print("here")

print(assemble(dot(u, n) * ds), assemble(dot(u, n) * ds(1)), assemble(dot(u, n) * ds(2)), assemble(dot(u, n) * ds(3)),
      assemble(dot(u, n) * ds(4)))

bcs[0].apply(u.vector())
print(assemble(dot(u, n) * ds), assemble(dot(u, n) * ds(1)), assemble(dot(u, n) * ds(2)), assemble(dot(u, n) * ds(3)),
      assemble(dot(u, n) * ds(4)))

bcs[1].apply(u.vector())
print(assemble(dot(u, n) * ds), assemble(dot(u, n) * ds(1)), assemble(dot(u, n) * ds(2)), assemble(dot(u, n) * ds(3)),
      assemble(dot(u, n) * ds(4)))

bcs[2].apply(u.vector())
print(assemble(dot(u, n) * ds), assemble(dot(u, n) * ds(1)), assemble(dot(u, n) * ds(2)), assemble(dot(u, n) * ds(3)),
      assemble(dot(u, n) * ds(4)))

bcs[3].apply(u.vector())
print(assemble(dot(u, n) * ds), assemble(dot(u, n) * ds(1)), assemble(dot(u, n) * ds(2)), assemble(dot(u, n) * ds(3)),
      assemble(dot(u, n) * ds(4)))

#bc.apply(u.vector())

#print(assemble(dot(u, n) * ds))

#print(u.compute_vertex_values())

