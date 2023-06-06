# -*- coding: utf-8 -*-


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os

set_log_level(LogLevel.WARNING)
#set_log_level(INFO)

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
G0 = Expression("1/(gam-1)*pow(1+alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])),gam)", alf=alf, bet=bet, gam=gam, degree=1)

# Define variational problem
V = FunctionSpace(mesh, 'CG', 1)
Q = FunctionSpace(mesh, "CG", 1)
#W = MixedFunctionSpace([P, P, P])
P = FiniteElement("P", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P, P, P]))

# Define unknown and test function(s)
(u_t, v_t, r_t) = TestFunctions(W)
w = Function(W)
(u, v, r) = (w[0], w[1], w[2])

# Solution at the previous time level
w_old = Function(W)
(u_old, v_old, r_old) = (w_old[0], w_old[1], w_old[2])

# Define boundary conditions
bcs = [DirichletBC(W.sub(0), Constant(0.0), boundaries, 1), \
       DirichletBC(W.sub(1), Constant(0.0), boundaries, 2)]

# Initial condition
u1 = project(Constant(0.0), V)
v1 = project(Constant(0.0), V)
r1 = project(R0, V)
assign(w.sub(0), u1)
assign(w.sub(1), v1)
assign(w.sub(2), r1)
w_old.assign(w)
# plot(r)
# interactive()

# Set the options for the time discretization
T = 5 # 5.
t = 0.
dt = 0.5
dti = 1.0 / dt
# dti = 2
# dt = 1 / dti

sm = r1 * dx()
sM = assemble(sm)
se = 0.5 * r1 * (u1 * u1 + v1 * v1) * dx() + 1. / (gam - 1) * r1 ** gam * dx()
sE = assemble(se)

tt = [t]
sMt = [sM]
sEt = [sE]

plt.figure(0)
dd = 0.00001
nN = 201
x = np.linspace(0., 5. - dd, nN)
gS = np.zeros(nN, 'float')
# for i in range(nN):
#    p = (x[i],0.)
#    gS[i] = r1(p)
# ss = "t = 0"
# plt.plot(x, gS, label=ss)

r00 = [r((0., 0.))]
rmin = [min(r1.vector())]
rmax = [max(r1.vector())]

# plot(r, key='r', title = " t=%g" % t)
File("result/full_" + str(0) + ".pvd") << r1

gT = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Define the variational formulation of the problem
F = dti * (r - r_old) * r_t * dx() + (r * u).dx(0) * r_t * dx() + (r * v).dx(1) * r_t * dx() \
    + dti * (r * u - r_old * u_old) * u_t * dx() + (r * u * u).dx(0) * u_t * dx() + (r * u * v).dx(
    1) * u_t * dx() + gam * r ** (gam - 1) * r.dx(0) * u_t * dx() \
    + dti * (r * v - r_old * v_old) * v_t * dx() + (r * v * u).dx(0) * v_t * dx() + (r * v * v).dx(
    1) * v_t * dx() + gam * r ** (gam - 1) * r.dx(1) * v_t * dx()

print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))

# u, v, r = w.split(True)
# plot(r)
# plt.show()

# Execute the time loop
while t < T - dt / 2:
    t += dt

    solve(F == 0, w, bcs, solver_parameters={"newton_solver":
                                                 {"relative_tolerance": 1e-15, "absolute_tolerance": 1e-15,
                                                  "relaxation_parameter": 1.0}})
    w_old.assign(w)
    u, v, r = w.split(True)
    # plot(r)
    # plt.show()

    sm = r * dx()
    sM = assemble(sm)
    se = 0.5 * r * (u * u + v * v) * dx() + 1. / (gam - 1) * r ** gam * dx()
    sE = assemble(se)

    tt.append(t)
    sMt.append(sM)
    sEt.append(sE)

    r00.append(r((0., 0.)))
    rmin.append(min(r.vector()))
    rmax.append(max(r.vector()))

    print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t, sM, sE))
    #    plot(r, key='r', title = " t=%g" % t)
    # for gt in gT:
    #     if t > gt - 0.1 * dt and t < gt + 0.1 * dt:
    #         for i in range(nN):
    #             p = (x[i], 0.)
    #             gS[i] = r(p)
    #         ss = "t = " + str(t)
    #         plt.plot(x, gS, label=ss)
    #         File("result/full_" + str(gt) + ".pvd") << r
    #         f1 = open("result/ref_" + str(gt) + ".npy", "wb")
    #         np.save(f1, gS)
    #         f1.close()

# f1 = open("result/ref_x.npy", "wb")
# np.save(f1, x)
# f1.close()
#
# f1 = open("fig2_3/tt.npy", "wb")
# np.save(f1, tt)
# f1.close()
# f1 = open("fig2_3/r00.npy", "wb")
# np.save(f1, r00)
# f1.close()
# f1 = open("fig2_3/rmin.npy", "wb")
# np.save(f1, rmin)
# f1.close()
# f1 = open("fig2_3/rmax.npy", "wb")
# np.save(f1, rmax)
# f1.close()
#
# plt.ylabel("$\\varrho$")
# plt.xlabel("$x_1$")
# plt.legend(loc=0)
# plt.grid()
#
# plt.figure(1)
# plt.plot(tt, r00, label="$r_0$")
# plt.plot(tt, rmin, label="$r_{min}$")
# plt.plot(tt, rmax, label="$r_{max}$")
# plt.ylabel("$r$")
# plt.xlabel("$t$")
# plt.legend(loc=0)
# plt.grid()
#
# plt.show()
# interactive()

