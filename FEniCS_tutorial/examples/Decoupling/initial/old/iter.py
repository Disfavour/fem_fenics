# -*- coding: utf-8 -*-


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os

set_log_level(CRITICAL)

#CRITICAL  = 50, // errors that may lead to data corruption and suchlike
#ERROR     = 40, // things that go boom
#WARNING   = 30, // things that may go boom later
#INFO      = 20, // information of general interest
#PROGRESS  = 16, // what's happening (broadly)
#TRACE     = 13, // what's happening (in detail)
#DBG       = 10  // sundry

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True

alf = 2.
bet = 20.0
gam = 1.4

# Mesh
nn = 200
mesh = RectangleMesh(Point(-5, -5),Point(5, 5), nn, nn, "right")
print mesh
#plot(mesh)
#interactive()

#isentropic vortex
boundaries = MeshFunction("uint", mesh, mesh.topology().dim() - 1)

class xBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and (x[0] < -5 + tol or x[0] > 5 - tol)
xB = xBoundary()
xB.mark(boundaries, 0)

class yBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and (x[1] < -5 + tol or x[1] > 5 - tol)
yB = yBoundary()
yB.mark(boundaries, 1)


ds = Measure("ds")[boundaries]

R0 = Expression("1+alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=1)

# Define variational problem
V = FunctionSpace(mesh, 'CG', 1)

# Define unknown and test function(s)
u = TrialFunction(V)
u_t = TestFunction(V)
v = TrialFunction(V)
v_t = TestFunction(V)
r = TrialFunction(V)
r_t = TestFunction(V)

# Solution at the previous time level
uk = Function(V)
vk = Function(V)
rk = Function(V)
r1 = Function(V)
u_old = Function(V)
v_old = Function(V)
r_old = Function(V)

# Define boundary conditions
bc1 = DirichletBC(V, Constant(0.0), boundaries, 0)
bc2 = DirichletBC(V, Constant(0.0), boundaries, 1)

# Initial condition
u_old = project(Constant(0.0), V)
v_old = project(Constant(0.0), V)
r_old = project(R0, V)

# Set the options for the time discretization
T = 5.
t = 0.
dt = 0.01
dti = 1.0/dt
dti = Constant(dti)
K = 1

sm = r_old * dx()
sM = assemble(sm)
se = 0.5*r_old*(u_old*u_old+v_old*v_old)*dx() + 1./(gam-1)*r_old**gam*dx()
sE = assemble(se)  

tt = [t]
sMt = [sM]
sEt = [sE]

plt.figure(0)
dd = 0.00001
nN = 201
x = np.linspace(0.,5.-dd, nN)
gS = np.zeros(nN, 'float')

r00 = [r_old((0.,0.))]
rmin = [min(r_old.vector())]
rmax = [max(r_old.vector())]

gT = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Define the variational formulation of the problem
F1 = dti*(r - r_old)*r_t*dx() + (r*uk).dx(0)*r_t*dx() + (r*vk).dx(1)*r_t*dx()
a1 = lhs(F1)
L1 = rhs(F1)

F2 = dti*(r1*u - r_old*u_old)*u_t*dx() + (r1*u*uk).dx(0)*u_t*dx() + (r1*u*vk).dx(1)*u_t*dx() + gam*r1**(gam-1)*r1.dx(0)*u_t*dx() 
a2 = lhs(F2)
L2 = rhs(F2)

F3 = dti*(r1*v - r_old*v_old)*v_t*dx() + (r1*v*uk).dx(0)*v_t*dx() + (r1*v*vk).dx(1)*v_t*dx() + gam*r1**(gam-1)*r1.dx(1)*v_t*dx()
a3 = lhs(F3)
L3 = rhs(F3)

# Execute the time loop
while t < T - dt/2:
    t += dt
    
    # Iteration
    rk.assign(r_old)
    uk.assign(u_old)
    vk.assign(v_old)
    for k in range(K):    
        solve(a1 == L1, r1)
        solve(a2 == L2, uk, bc1)
        solve(a3 == L3, vk, bc2)
        rk.assign(r1)
    r_old.assign(rk)
    u_old.assign(uk)
    v_old.assign(vk)
    
    sm = r_old * dx()
    sM = assemble(sm)
    se = 0.5*r_old*(u_old*u_old+v_old*v_old)*dx() + 1./(gam-1)*r_old**gam*dx()
    sE = assemble(se)   
    
    tt.append(t)
    sMt.append(sM)
    sEt.append(sE)

    r00.append(r_old((0.,0.)))
    rmin.append(min(r_old.vector()))
    rmax.append(max(r_old.vector()))
    
    print("conservation: t %0.3f\t M %0.10f\t E %0.10f" % (t,sM,sE))
#    plot(r_old, key='r', title = " t=%g" % t)
    for gt in gT:
        if t > gt - 0.1*dt and t < gt + 0.1*dt:    
            for i in range(nN):
                p = (x[i],0.)
                gS[i] = r_old(p)
            ss = "t = " + str(t)  
            plt.plot(x, gS, label=ss)
            f1 = open("K_" + str(K) + "_" + str(gt) + ".npy", "wb")
            np.save(f1, gS)
            f1.close()
            
f1 = open("K_" + str(K) + "_" + str(dt) + "_t.npy", "wb")
np.save(f1, tt)
f1.close()
f1 = open("K_" + str(K) + "_" + str(dt) + "_E.npy", "wb")
np.save(f1, sEt)
f1.close()
            
plt.ylabel("$\\varrho$")
plt.xlabel("$x_1$")  
plt.legend(loc=0)
plt.grid()

plt.show()

