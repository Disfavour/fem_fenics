from fenics import *
import numpy as np
from scipy.constants import g
from scipy.integrate import quad
import matplotlib.pyplot as plt

set_log_level(LogLevel.WARNING)


hl, hr = 10, 1
h1 = 3.9618#4.0612#3.9618
u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
D1 = - (g * hl) ** 0.5
D2 = u1 - (g*h1) ** 0.5
D3 = u1 * h1 / (h1 - hr)
domain_size = 5


def solution_h(x, t):
    'x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))'
    if x < D1*t:
        return hl
    elif x < D2*t:
        return 1/(9*g) * (2*(g*hl) ** 0.5 - x/t) ** 2
    elif x < D3*t:
        return h1
    else:
        return hr

def solution_u(x, t):
    'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'
    if x < D1*t:
        return 0
    elif x < D2*t:
        return 1/3 * (2*(g*hl)**0.5 + 2*x/t)
    elif x < D3*t:
        return u1
    else:
        return 0


def calc_E(x, t):
    h = solution_h(x, t)
    u = solution_u(x, t)
    return 0.5 * (h * u ** 2 + g * h ** 2)


def analytic(ts):
    E = []
    for t in ts:
        E.append(quad(calc_E, -5, 5, (t,), limit=200)[0])
    return np.array(E)


def sw_1d(mesh_size=200, tau=0.005):
    degree = 1
    T = 0.5

    ms, Es = [], []

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    U = FiniteElement('P', mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    w = Function(W)
    h, u = split(w)

    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    def collect_data():
        ms.append(assemble(m_eq))
        Es.append(assemble(E_eq))
        print(f'Time {t:>7.5f} m {ms[-1]:>7.5f} E {Es[-1]:>7.5f}')

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=1), W))
    w.assign(project(w_exact, W))
    collect_data()

    for t in ts[1:]:
        w_exact.t = t
        w.assign(project(w_exact, W))
        collect_data()

    return map(np.array, (ts, Es))


def plot_E(ts, Es, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    for E in Es:
        plt.plot(ts, E)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0], ts[-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()


if __name__ == '__main__':
    mesh_size = 400
    tau = 0.005

    ts, Es = sw_1d(mesh_size, tau)
    E_analytic = analytic(ts)

    plot_E(ts, (Es, E_analytic), ('fem', 'analytic'), f'E1_ms{mesh_size}_tau{tau}.png')
