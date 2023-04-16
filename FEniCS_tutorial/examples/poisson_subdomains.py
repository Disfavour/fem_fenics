from fenics import *
import numpy as np
import matplotlib.pyplot as plt


tol = 1E-14
k_0 = 1.0
k_1 = 0.01


mesh = UnitSquareMesh(20, 20)

V = FunctionSpace(mesh, "P", 1)

materials = MeshFunction("size_t", mesh, mesh.topology().dim())

subdomain_0 = CompiledSubDomain("x[1] <= 0.5 + tol", tol=tol)
subdomain_1 = CompiledSubDomain("x[1] >= 0.5 - tol", tol=tol)

materials.set_all(0)
subdomain_1.mark(materials, 1)

#File("materials.xml") << materials

# plot(materials)
# plt.show()

# subdomain_3 = CompiledSubDomain("pow((x[0] - 0.5), 2) + pow((x[1] - 0.5), 2) <= 0.2", tol=tol)
#
# cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
# cell_markers.set_all(False)
# subdomain_3.mark(cell_markers, True)
#
# mesh = refine(mesh, cell_markers)
#
# plot(mesh)
# plt.show()


class K(UserExpression):
    def __init__(self, materials, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.materials = materials
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.materials[cell.index] == 0:
            values[0] = self.k_0
        else:
            values[0] = self.k_1


kappa = K(materials, k_0, k_1, degree=0)

# boundary_R = CompiledSubDomain("on_boundary && near(x[0], 1, tol)", tol=1E-14)

bcs = DirichletBC(V, Constant(0), "on_boundary")


# exact solution
u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = kappa*dot(grad(u), grad(v))*dx
L = f*v*dx

u = Function(V)
solve(a == L, u, bcs)

plt.colorbar(plot(u))
plot(mesh)

# vtkfile = File("poisson/solution.pvd")
# vtkfile << u

error_L2 = errornorm(u_D, u, "L2")

vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

print("error_L2  =", error_L2)
print("error_max =", error_max)

plt.show()
