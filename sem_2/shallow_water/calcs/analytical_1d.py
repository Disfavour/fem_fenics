from fenics import *
import numpy as np

# Создание сетки
mesh = UnitSquareMesh(32, 32)

# Определение функционального пространства
V = FunctionSpace(mesh, "Lagrange", 1)

# Определение функций и тестовых функций
u = Function(V)
v = TestFunction(V)

# Определение параметров уравнения
k = Constant(1.0)
f = Expression("x[0] * sin(x[1])", degree=2)

# Определение нелинейной формы F
F = inner((1 + u**2)*grad(u), grad(v)) * dx - f * v * dx

# Определение граничных условий
u_D = Expression("x[0] * x[0] + x[1] * x[1]", degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Определение нелинейного решателя методом Ньютона
problem = NonlinearVariationalProblem(F, u, bc, J=derivative(F, u))
solver = NonlinearVariationalSolver(problem)

# Параметры решателя Ньютона
# solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-6
# solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
# solver.parameters["newton_solver"]["maximum_iterations"] = 25
# solver.parameters["newton_solver"]['relaxation_parameter'] = 1.0

# #print(list_linear_solver_methods())
# iterative_solver = True

# if iterative_solver:
#     solver.parameters['nonlinear_solver']='newton'
#     #solver.parameters['newton_solver']['linear_solver'] = 'mumps' 'gmres'

#     solver.parameters['newton_solver']['linear_solver'] = 'gmres'
#     solver.parameters['newton_solver']['preconditioner'] = 'ilu'
#     solver.parameters['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-9
#     solver.parameters['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-7
#     solver.parameters['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
#     #solver.parameters['newton_solver']['krylov_solver']['gmres']['restart'] = 40
#     #solver.parameters['newton_solver']['krylov_solver']['preconditioner']['ilu']['fill_level'] = 0

#     # solver.parameters['linear_solver'] = 'gmres'
#     # solver.parameters['preconditioner'] = 'ilu'
#     # solver.parameters['krylov_solver']['absolute_tolerance'] = 1E-9
#     # solver.parameters['krylov_solver']['relative_tolerance'] = 1E-7
#     # solver.parameters['krylov_solver']['maximum_iterations'] = 1000
#     # solver.parameters['krylov_solver']['gmres']['restart'] = 40
#     # solver.parameters['krylov_solver']['preconditioner']['ilu']['fill_level'] = 0

# #set_log_level(LogLevel.PROGRESS)

J = derivative(F, u)

problem = NonlinearVariationalProblem(F, u , bc, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['nonlinear_solver'] = 'snes'
prm['snes_solver']['line_search'] = 'bt'
prm['snes_solver']['linear_solver'] = 'mumps'
#prm[“snes_solver”][“linear_solver”] = “gmres”
#prm[“snes_solver”][“preconditioner”] = “hypre_amg”
#prm[“snes_solver”][“krylov_solver”][“nonzero_initial_guess”] = False
prm['snes_solver']['report'] = True

# Решение нелинейной задачи
solver.solve()