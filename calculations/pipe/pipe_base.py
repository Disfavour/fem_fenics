from fenics import *
import numpy as np
from scipy.constants import g, pi, R
import matplotlib.pyplot as plt


set_log_level(LogLevel.WARNING)

P_in, P_out, m_in, m_out = [], [], [], []

t_max = 24 * 3600
mesh_size = 400
tau = 0.1

L = 1e5
D = 0.6
T = 278
h = 0

eps = 0.000617
Re = 2000
A = pi * D**2 / 4
f = (-2*np.log10(eps/D/3.7 - 4.518/Re*np.log10(6.9/Re + (eps/D/3.7)**1.11))) ** -2 

S = 0.6
M_air = 28.964917 / 1000
M = S * M_air
Rs = R / M

#k = 0
Z = 1

sigma = 1

ts = np.arange(0, t_max+tau/2, tau)
mesh = IntervalMesh(mesh_size, 0, L)

RHO = FiniteElement('P', mesh.ufl_cell(), 1)
U = FiniteElement('P', mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([RHO, U]))

u_out = Expression('(t < 2 ? 20*t+70 : (t < 10 ? 110 : (t < 12 ? -40*t+510 : (t < 22 ? 30 : 20*t-410))))/A', t=0, A=A, degree=1)

bc = [DirichletBC(W.sub(0), Constant(5e6 / (Z*Rs*T)), 'on_boundary && near(x[0], 0)'),
      DirichletBC(W.sub(1), u_out, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

w, wn = Function(W), Function(W)
rt, ut = TestFunctions(W)
r, u = split(w)
rn, un = split(wn)

r = np.array((r, rn))
u = np.array((u, un))

dt = lambda f: (f[0]-f[1])/tau
s = lambda f: sigma*f[0] + (1-sigma)*f[1]

F = (dt(r) + s(r*u).dx(0)) * rt*dx \
    + (dt(r*u) + s(r*u**2 + r*Z*Rs*T).dx(0) + s(r*u*abs(u))*f/(2*D) + s(r)*g*h/L) * ut*dx

def collect_data():
    P_in.append(w.sub(0)(0) * Z*Rs*T)
    P_out.append(w.sub(0)(L) * Z*Rs*T)
    # Q
    m_in.append(w.sub(1)(0) * A)
    m_out.append(w.sub(1)(L) * A)

    print(f'Time {t:>7.5f} {w.sub(0)(0):>7.5f} {w.sub(0)(L):>7.5f} {P_in[-1]/1e5:>7.5f} {P_out[-1]/1e5:>7.5f} {m_in[-1]} {m_out[-1]}')

    if t % 600 == 0:
        print(f'Time {t:>7.5f}')

        plt.figure()
        plt.plot(ts[:len(P_in)], P_in)
        plt.savefig("P_in")
        plt.close()

        plt.figure()
        plt.plot(ts[:len(P_out)], P_out)
        plt.savefig("P_out")
        plt.close()

        plt.figure()
        plt.plot(ts[:len(m_in)], m_in)
        plt.savefig("m_in")
        plt.close()

        plt.figure()
        plt.plot(ts[:len(m_out)], m_out)
        plt.savefig("m_out")
        plt.close()

t = 0
w.assign(project(Expression(('(5e6-8*x[0])/(Z*Rs*T)', '70/A'), Z=Z, Rs=Rs, T=T, A=A, degree=1), W))
# (5e6-8*x[0])/(Z*Rs*T)
collect_data()
wn.assign(w)

# plot(w.sub(0))
# plt.figure()
# plot(w.sub(1))
# plt.show()
# exit()

for t in ts[1:]:
    u_out.t = t / 3600

    solve(F == 0, w, bc)

    collect_data()

    wn.assign(w)
