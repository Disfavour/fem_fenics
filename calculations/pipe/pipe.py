from fenics import *
import numpy as np
from scipy.constants import g, pi, R
import matplotlib.pyplot as plt


set_log_level(LogLevel.WARNING)

P_in, P_out, m_in, m_out = [], [], [], []

t_max = 24 * 3600
L = 1e5

k = 0
S = 0.6
D = 0.6
T = 278
h = 0
eps = 0.000617
Re = 2000

A = pi * D**2 / 4
f = (-2*np.log10(eps/D/3.7 - 4.518/Re*np.log10(6.9/Re + (eps/D/3.7)**1.11))) ** -2 



M_air = 28.964917 / 1000
M = S * M_air
Rs = R / M

rho = 4.2e6 / (Rs * T) #rho = 0.73

sigma = 1

mesh_size = 5
tau = 0.1

ts = np.arange(0, t_max+tau/2, tau)
mesh = IntervalMesh(mesh_size, 0, L)

P_ = FiniteElement('P', mesh.ufl_cell(), 1)
m_ = FiniteElement('P', mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P_, m_]))

m_out_expr = Expression('rho*(t < 2 ? 20*t+70 : (t < 10 ? 110 : (t < 12 ? -40*t+510 : (t < 22 ? 30 : 20*t-410))))', t=0, rho=rho, degree=1)

bc = [DirichletBC(W.sub(0), Constant(5e6), 'on_boundary && near(x[0], 0)'),
      DirichletBC(W.sub(1), m_out_expr, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

w, wn = Function(W), Function(W)
Pt, mt = TestFunctions(W)
P, m = split(w)
Pn, mn = split(wn)

P = np.array((P, Pn))
m = np.array((m, mn))

dt = lambda f: (f[0]-f[1])/tau
s = lambda f: sigma*f[0] + (1-sigma)*f[1]

Z = 1 #+ k*s(P)

F = (dt(P/Z)/(Rs*T) + s(m).dx(0)/A) * Pt*dx \
    + (dt(m)/A + Rs*T/A**2*s(Z*m**2/P).dx(0) + s(P).dx(0) + Rs*T*f/(2*D*A**2)*s(m*abs(m)*Z/P) + g*h/(L*Rs*T)*s(P/Z)) * mt*dx

def collect_data():
    P_in.append(w.sub(0)(0) / 1e6)
    P_out.append(w.sub(0)(L) / 1e6)
    m_in.append(w.sub(1)(0) / rho)
    m_out.append(w.sub(1)(L) / rho)

    print(P_out[-1], rho)

    if True:#t % 600 == 0:
        print(f'Time {t:>7.5f} rho {rho:>7.5f}')

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
w.assign(project(Expression(('1e5*42', 'rho*70'), rho=rho, degree=1), W))
# 1e5*(50-0.0008*x[0]) 1e5*42
collect_data()
wn.assign(w)



# plot(w.sub(0))
# plt.figure()
# plot(w.sub(1))
# plt.show()
# exit()



for t in ts[1:]:
    m_out_expr.t = t / 3600
    m_out_expr.rho = rho

    solve(F == 0, w, bc)

    #rho = w.sub(0).compute_vertex_values().mean() / (Z*Rs*T)

    collect_data()

    wn.assign(w)

plt.plot(ts, P_in)

plt.figure()
plt.plot(ts, P_out)

plt.figure()
plt.plot(ts, m_in)

plt.figure()
plt.plot(ts, m_out)

plt.show()
