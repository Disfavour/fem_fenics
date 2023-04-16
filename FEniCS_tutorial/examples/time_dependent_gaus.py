from fenics import *
import matplotlib.pyplot as plt


T = 2.0
num_steps = 50
dt = T / num_steps

nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)

V = FunctionSpace(mesh, "P", 1)

bc = DirichletBC(V, Constant(0), "on_boundary")

u_0 = Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))", degree=2, a=5)
u_n = interpolate(u_0, V)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

u = Function(V)
t = 0
for n in range(num_steps):
    # Обновляем время
    t += dt
    # Считаем u на следующем шаге
    solve(a == L, u, bc)

    plot(u)
    plt.show()

    # Обновляем результат предыдущего шага
    u_n.assign(u)

