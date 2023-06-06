# -*- coding: utf-8 -*-


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os

from ufl import nabla_div

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
#parameters["form_compiler"]["cpp_optimize"] = True

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
V = VectorFunctionSpace(mesh, "P", 1)
S = FunctionSpace(mesh, 'CG', 1)

# u_ = TrialFunction(V)
# v_ = TrialFunction(V)
u__ = TrialFunction(V)
u_, v_ = split(u__)
r_ = TrialFunction(S)

# ut = TestFunction(V)
# vt = TestFunction(V)
ut_ = TestFunction(V)
ut, vt = split(ut_)
rt = TestFunction(S)

# u = Function(V)
# v = Function(V)
ur = Function(V)
u, v = split(ur)
r = Function(S)

# uk = Function(V)
# vk = Function(V)
uk_ = Function(V)
uk, vk = split(uk_)
rk = Function(S)

# un = Function(V)
# vn = Function(V)
un_ = Function(V)
un, vn = split(un_)
rn = Function(S)

# Define boundary conditions
bc1 = DirichletBC(V.sub(0), Constant(0.0), boundaries, 1)
bc2 = DirichletBC(V.sub(1), Constant(0.0), boundaries, 2)

# Initial condition
# un = project(Constant(0.0), V)
# vn = project(Constant(0.0), V)
rn = project(R0, S)

# Set the options for the time discretization
T = 1  # 5.
t = 0.
dt = 0.02
dti = 1.0 / dt
# dti = 35
# dt = 1 / dti
dti = Constant(dti)
K = 1

sm = rn * dx()
se = 0.5 * rn * (un * un + vn * vn) * dx() + 1. / (gam - 1) * rn ** gam * dx()
sM = assemble(sm)
sE = assemble(se)

# Define the variational formulation of the problem
F1 = dti * (r_ - rn) * rt * dx() + (r_ * uk).dx(0) * rt * dx() + (r_ * vk).dx(1) * rt * dx()
a1 = lhs(F1)
L1 = rhs(F1)

F2 = dti * (r * u_ - rn * un) * ut * dx() + (r * u_ * uk).dx(0) * ut * dx() + (r * u_ * vk).dx(
    1) * ut * dx() + gam * r ** (gam - 1) * r.dx(0) * ut * dx()
a2 = lhs(F2)
L2 = rhs(F2)

F3 = dti * (r * v_ - rn * vn) * vt * dx() + (r * v_ * uk).dx(0) * vt * dx() + (r * v_ * vk).dx(
    1) * vt * dx() + gam * r ** (gam - 1) * r.dx(1) * vt * dx()
a3 = lhs(F3)
L3 = rhs(F3)

F23 = dot((r*u__-rn*un_)*dti, ut_) * dx\
    + dot(nabla_div(outer(uk_, r*u__)), ut_) * dx\
    + dot(grad(r**gam), ut_) * dx
#F2 + F3
a23 = lhs(F23)
L23 = rhs(F23)

print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))

uk_ = Function(V)

while t < T - dt / 2:
    t += dt

    rk.assign(rn)
    # uk.assign(un)
    # vk.assign(vn)
    uk_.assign(un_)

    for k in range(K):
        solve(a1 == L1, r)
        # solve(a2 == L2, u, bc1)
        # solve(a3 == L3, v, bc2)
        solve(a23 == L23, ur, [bc1, bc2])

        # uk.assign(u)
        # vk.assign(v)
        uk_.assign(ur)
        rk.assign(r)

    rn.assign(rk)
    # un.assign(uk)
    # vn.assign(vk)
    un_.assign(uk_)

    sM = assemble(sm)
    sE = assemble(se)
    print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))
