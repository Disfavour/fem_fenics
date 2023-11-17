import matplotlib.pyplot as plt
from dolfin import *
import numpy as np


class Problem(NonlinearProblem):
    def __init__(self, J, F, bcs):
        self.bilinear_form = J
        self.linear_form = F
        self.bcs = bcs
        NonlinearProblem.__init__(self)

    def F(self, b, x):
        assemble(self.linear_form, tensor=b)
        for bc in self.bcs:
            bc.apply(b, x)

    def J(self, A, x):
        assemble(self.bilinear_form, tensor=A)
        for bc in self.bcs:
            bc.apply(A)


class CustomSolver(NewtonSolver):
    def __init__(self):
        NewtonSolver.__init__(self)#,# mesh.mpi_comm(),
                              #PETScKrylovSolver(), PETScFactory.instance())

    def solver_setup(self, A, P, problem, iteration):
        self.linear_solver().set_operator(A)

        PETScOptions.set("ksp_type", "gmres")
        PETScOptions.set("ksp_monitor")
        PETScOptions.set("pc_type", "ilu")

        self.linear_solver().set_from_options()


mesh = UnitSquareMesh(32, 32)

V = FunctionSpace(mesh, "CG", 1)
g = Constant(1.0)
bcs = [DirichletBC(V, g, "near(x[0], 1.0) and on_boundary")]
u = Function(V)
v = TestFunction(V)
f = Expression("x[0]*sin(x[1])", degree=2)
F = inner((1 + u**2)*grad(u), grad(v))*dx - f*v*dx
J = derivative(F, u)
u_new = Function(V)
u_new.assign(u)

# problem = Problem(J, F, bcs)
# custom_solver = CustomSolver()
# custom_solver.solve(problem, u.vector())

J = derivative(F, u)

problem = NonlinearVariationalProblem(F, u , bcs, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['nonlinear_solver'] = 'snes'
prm['snes_solver']['line_search'] = 'bt'
prm['snes_solver']['linear_solver'] = 'mumps'
#prm[“snes_solver”][“linear_solver”] = “gmres”
#prm[“snes_solver”][“preconditioner”] = “hypre_amg”
#prm[“snes_solver”][“krylov_solver”][“nonzero_initial_guess”] = False
prm['snes_solver']['report'] = True

solver.solve()


plt.figure()
plot(u, title="Solution")

plt.figure()
plot(u_new)

plt.figure()
plot(grad(u), title="Solution gradient")

plt.show()