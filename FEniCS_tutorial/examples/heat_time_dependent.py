from fenics import *
import numpy as np
import matplotlib.pyplot as plt


T = 2.0 # final time
num_steps = 10 # number of time steps
dt = T / num_steps # time step size
alpha = 3 # parameter alpha
beta = 1.2 # parameter beta

# Create mesh and define function space
nx = ny = 8
mesh = UnitSquareMesh(nx, ny)

V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
u_D = Expression("1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t", degree=2, alpha=alpha, beta=beta, t=0)

bc = DirichletBC(V, u_D, "on_boundary")

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*alpha)
F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):
    # Update current time
    t += dt
    u_D.t = t
    # Compute solution
    solve(a == L, u, bc)
    # Plot solution
    print(errornorm(u_D, u))
    plot(u)
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
    print(error_max)
    #plt.show()
    # Update previous solution
    u_n.assign(u)

