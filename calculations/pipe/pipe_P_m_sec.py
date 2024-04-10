from fenics import *
import numpy as np
from scipy.constants import pi, R
import matplotlib.pyplot as plt


set_log_level(LogLevel.WARNING)

P_in, P_out, m_in, m_out = [], [], [], []

t_max = 24 * 3600
mesh_size = 500
tau = 0.1

L = 1e5
T = 278

D = 0.6
A = pi * D**2 / 4

eps = 0.000617
Re = 2000
f = (-2*np.log10(eps/D/3.7 - 4.518/Re*np.log10(6.9/Re + (eps/D/3.7)**1.11))) ** -2 

S = 0.6
M_air = 28.964917 / 1000
M = S * M_air
Rs = R / M

# sec
P_left = 5e6
P_0_right = 4.2e6
Q_0 = 70
# t в часах
m_out_expr = Expression('rho_out * (t < 2 ? 20*t+70 : (t < 10 ? 110 : (t < 12 ? -40*t+510 : (t < 22 ? 30 : 20*t-410))))', rho_out=0, t=0, degree=1)

ts = np.arange(0, t_max+tau/2, tau)
mesh = IntervalMesh(mesh_size, 0, L)

P_ = FiniteElement('P', mesh.ufl_cell(), 1)
m_ = FiniteElement('P', mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P_, m_]))

bc = [DirichletBC(W.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)'),
      DirichletBC(W.sub(1), m_out_expr, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

w, wn = Function(W), Function(W)
Pt, mt = TestFunctions(W)
P, m = split(w)
Pn, mn = split(wn)

F = ((P-Pn)/(tau*Rs*T) + m.dx(0)/A) * Pt*dx \
    + ((m-mn)/(tau*A) + Rs*T/A**2*(m**2/P).dx(0) + P.dx(0) + Rs*T*f*m*abs(m)/(2*D*A**2*P)) * mt*dx

def collect_data():
    P_in.append(w.sub(0)(0) / 1e5)
    P_out.append(w.sub(0)(L) / 1e5)
    m_in.append(w.sub(1)(0))
    m_out.append(w.sub(1)(L))

    if t % 100 == 0:
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
w.assign(project(Expression(('(P_left + x[0]*(P_right - P_left)/L)', '(P_left + x[0]*(P_right - P_left)/L) / (Rs*T) * Q_0'),
                            P_left=P_left, P_right=P_0_right, L=L, Rs=Rs, T=T, Q_0=Q_0, degree=1), W))
# w.assign(project(Expression(('(P_left + x[0]*(P_right - P_left)/L)', '0.73 * Q_0'),
#                             P_left=P_left, P_right=P_0_right, L=L, Rs=Rs, T=T, Q_0=Q_0, degree=1), W))
collect_data()
wn.assign(w)



# plot(w.sub(0))
# plt.figure()
# plot(w.sub(1) / (w.sub(0)(L)/(Rs*T)))
# plt.show()
# exit()



for t in ts[1:]:
    m_out_expr.t = t / 3600
    m_out_expr.rho_out = w.sub(0)(L)/(Rs*T)
    print(m_out_expr.rho_out)
    

    solve(F == 0, w, bc)

    collect_data()

    wn.assign(w)
