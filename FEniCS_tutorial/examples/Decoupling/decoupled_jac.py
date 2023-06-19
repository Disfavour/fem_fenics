# -*- coding: utf-8 -*-


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os

set_log_level(LogLevel.WARNING)
#set_log_level(CRITICAL)

# CRITICAL  = 50, // errors that may lead to data corruption and suchlike
# ERROR     = 40, // things that go boom
# WARNING   = 30, // things that may go boom later
# INFO      = 20, // information of general interest
# PROGRESS  = 16, // what's happening (broadly)
# TRACE     = 13, // what's happening (in detail)
# DBG       = 10  // sundry

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

alf = 2.
bet = 20.0
gam = 1.4

# Mesh
nn = 200
mesh = RectangleMesh(Point(-5, -5), Point(5, 5), nn, nn, "right")
# print
# mesh
# plot(mesh)
# interactive()

# isentropic vortex
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)


class xBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and (x[0] < -5 + tol or x[0] > 5 - tol)


xB = xBoundary()
xB.mark(boundaries, 1)


class yBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and (x[1] < -5 + tol or x[1] > 5 - tol)


yB = yBoundary()
yB.mark(boundaries, 2)

#ds = Measure("ds")[boundaries]

R0 = Expression("1+alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=1)

# Define variational problem
V = FunctionSpace(mesh, 'CG', 1)

# Define unknown and test function(s)
ut, vt, rt = TestFunction(V), TestFunction(V), TestFunction(V)
u_, v_, r_ = TrialFunction(V), TrialFunction(V), TrialFunction(V)

# Solution at the previous time level
uk, vk, rk = Function(V), Function(V), Function(V)
un, vn, rn = Function(V), Function(V), Function(V)
r = Function(V)

# Define boundary conditions
bc1 = DirichletBC(V, Constant(0.0), boundaries, 1)
bc2 = DirichletBC(V, Constant(0.0), boundaries, 2)

# Initial condition
un = project(Constant(0.0), V)
vn = project(Constant(0.0), V)
rn = project(R0, V)

# Set the options for the time discretization
T = 1  # 5.
t = 0.
dt = 0.02
dti = 1.0 / dt
tau = dt
# dti = 35
# dt = 1 / dti
dti = Constant(dti)
K = 1

Es = []

sm = rn * dx()
sM = assemble(sm)
se = 0.5 * rn * (un * un + vn * vn) * dx() + 1. / (gam - 1) * rn ** gam * dx()
sE = assemble(se)

Es.append(sE)

# Define the variational formulation of the problem
F1 = (r_-rn)/tau * rt * dx\
    + (r_*uk).dx(0) * rt * dx + (r_*vk).dx(1) * rt * dx

a1, L1 = lhs(F1), rhs(F1)

F2 = (r*u_ - rn*un)/tau * ut * dx\
     + (r*u_*uk).dx(0) * ut * dx + (r*u_*vk).dx(1) * ut * dx\
     + gam * r ** (gam - 1) * r.dx(0) * ut * dx
a2, L2 = lhs(F2), rhs(F2)

F3 = (r*v_ - rn*vn)/tau * vt * dx\
     + (r*v_*uk).dx(0) * vt * dx + (r*v_*vk).dx(1) * vt * dx\
     + gam * r ** (gam - 1) * r.dx(1) * vt * dx
a3, L3 = lhs(F3), rhs(F3)

print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))

# u, v, r = w.split(True)
# plot(r)
# plt.show()

uk_ = Function(V)

# Execute the time loop
while t < T - dt / 2:
    t += dt

    # Iteration
    rk.assign(rn)
    uk.assign(un)
    vk.assign(vn)
    for k in range(K):
        solve(a1 == L1, r)
        solve(a2 == L2, uk_, bc1)
        solve(a3 == L3, vk, bc2)
        rk.assign(r)
        uk.assign(uk_)
    rn.assign(rk)
    un.assign(uk)
    vn.assign(vk)

    sm = rn * dx()
    sM = assemble(sm)
    se = 0.5 * rn * (un * un + vn * vn) * dx() + 1. / (gam - 1) * rn ** gam * dx()
    sE = assemble(se)

    Es.append(sE)

    print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))

np.save("decoupled_jac2.npy", np.array(Es))
