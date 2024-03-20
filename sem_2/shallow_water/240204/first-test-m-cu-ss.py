from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# set_log_level(CRITICAL)

# first order

test = 1
m = 200
cu = 0.25
ss = 0.125

gam = 0.5


# Mesh
mesh = IntervalMesh(m, 0., 1.)

p = 3
V = FunctionSpace(mesh, 'CG', p)

l1 = 0.05
l2 = 0.25
if test == 1:
    U0 = Expression("l1 < x[0] && x[0] < l2? 1 : 0", l1=l1, l2=l2, degree=p+1)

# Define unknown and test function(s)
u = Function(V)
u_t = TestFunction(V)

# Solution at the previous time level
u_old = Function(V)

# # Initial condition
# u0 = project(U0, V)
u0 = interpolate(U0, V)
plot(u0)
plt.show()
u_old.assign(u0)
u.assign(u_old)


# time discretization
T = 0.75
t = 0.
dt0 = 1./m
dt = cu*dt0

# local
# ss = 0.01*abs(u.dx(0))


# Define the variational formulation of the problem
F = 1.0/dt*(u - u_old)*u_t*dx + 0.5*(u.dx(0)+u_old.dx(0))*u_t*dx # CN
# F = F + ss*(u.dx(0)-u_old.dx(0))*u_t*dx # схема с весом ss = sig - 0.5
F = F + ss*dt*u.dx(0)*u_t.dx(0)*dx  # искуственная вязкость  

bc = DirichletBC(V, Constant(0.), 'on_boundary && near(x[0],0)')

# Execute the time loop
while t < T - dt/2:
    t += dt

    solve(F == 0, u, bc)
    u_old.assign(u)
            
    if abs(t - 0.25) < 0.1*dt:
        plot(u)
        plt.show()
    
    if abs(t - 0.5) < 0.1*dt:
        plot(u)
        plt.show()
            
    if abs(t - 0.75) < 0.1*dt:
        plot(u)
        plt.show()
        






