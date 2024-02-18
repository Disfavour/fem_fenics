from fenics import *
import numpy as np
from scipy.constants import g
import scipy
import matplotlib.pyplot as plt
from os.path import dirname, join


dpi=300

def plot_Es(ts, Es, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for t, E in zip(ts, Es):
        plt.plot(t, E)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_Es_dif(ts, Es, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for t, E in zip(ts[:-1], Es[:-1]):
        plt.plot(t, E - Es[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E_{fem} - E_{quad}$')
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_Es_tau(ts, Es, names, ylabel, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for t, E in zip(ts, Es):
        plt.plot(t, E)
    plt.xlabel(r'$t$')
    plt.ylabel(ylabel)
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_m(t, me, m_sc, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, me)
    plt.plot(t, m_sc)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$m$')
    plt.xlim(t[0], t[-1])
    plt.legend(('fem', 'quad'))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_E(t, Ee, E_sc, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, Ee)
    plt.plot(t, E_sc)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(t[0], t[-1])
    plt.legend(('fem', 'quad'))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_E1(t, E, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, E)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\frac {1} {2} h u^2$')
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_E2(t, E, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, E)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\frac {1} {2} g h^2$')
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()


def plot_h(x, hs, legend, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for h in hs:
        plt.plot(x, h)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$h$')
    plt.xlim(x[0], x[-1])
    plt.grid()
    plt.legend(legend)
    plt.savefig(fname)
    plt.close()


set_log_level(LogLevel.WARNING)

hl, hr = 10, 1
h1 = 3.9618
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

def calc_m(x, t):
    return solution_h(x, t)

def calc_E(x, t):
    h = solution_h(x, t)
    u = solution_u(x, t)
    return 0.5 * (h * u ** 2 + g * h ** 2)

def calc_E1(x, t):
    h = solution_h(x, t)
    u = solution_u(x, t)
    return 0.5 * h * u ** 2

def calc_E2(x, t):
    h = solution_h(x, t)
    return 0.5 * g * h ** 2

def check(ts):
    # m_eq = h * dx
    # E_eq = 0.5*(h*u*u + g*h*h) * dx
    m, E = [], []
    E1, E2 = [], []
    for t in ts:
        m.append(scipy.integrate.quad(calc_m, -5, 5, (t,), limit=200)[0])
        E.append(scipy.integrate.quad(calc_E, -5, 5, (t,), limit=200)[0])
        E1.append(scipy.integrate.quad(calc_E1, -5, 5, (t,), limit=200)[0])
        E2.append(scipy.integrate.quad(calc_E2, -5, 5, (t,), limit=200)[0])
    return m, E, E1, E2

def sw_1d(mesh_size=200, tau=0.005, s=1, time_moments=[0.0, 0.25, 0.5]):
    degree = 1
    T = 0.5

    m, E, me, Ee = [], [], [], []
    res_h = []

    # hl, hr = 10, 1
    
    # h1 = 3.9618
    # u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    # D1 = - (g * hl) ** 0.5
    # D2 = u1 - (g*h1) ** 0.5
    # D3 = u1 * h1 / (h1 - hr)
    # domain_size = 5

    t = 0
    time_steps = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    x = mesh.coordinates().flatten()

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    we, wne = Function(W), Function(W)
    he, ue = split(we)
    hne, une = split(wne)

    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    hse = s*he + (1-s)*hne
    use = s*ue + (1-s)*une

    # F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * ut*dx
    
    F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0)) * ut*dx - g/2*(hs*hs) * ut.dx(0) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    m_eq_exact = he * dx
    E_eq_exact = 0.5*(he*ue*ue + g*he*he) * dx

    def collect_data():
        #m, E, me, Ee = map(assemble, (m_eq, E_eq, m_eq_exact, E_eq_exact))
        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        me.append(assemble(m_eq_exact))
        Ee.append(assemble(E_eq_exact))
        print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

        if np.isclose(t, time_moments).any():
            res_h.append(w.sub(0).compute_vertex_values())

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=1), W))

    w_exact.t = t
    we.assign(project(w_exact, W))

    collect_data()

    wn.assign(w)
    wne.assign(we)

    for t in time_steps[1:]:
        solve(F == 0, w, bc, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-8,
            'relative_tolerance': 1e-8,
            'maximum_iterations': 50,
            'relaxation_parameter': 1.0,
        }})

        w_exact.t = t
        we.assign(project(w_exact, W))

        collect_data()

        wn.assign(w)
        wne.assign(we)

    return map(np.array, (time_steps, m, E, me, Ee, x, res_h))


if __name__ == '__main__':
    base_dir = dirname(__file__)
    images_dir = join(base_dir, 'images')

    ms=400
    tau=0.00125
    s=1

    time_moments=[0.0, 0.25, 0.5]
    ts, m, Es, me, Ee, x, hs = sw_1d(ms, tau, s, time_moments)

    plot_h(x, hs, ('$t=0$', '$t=0.25$', '$t=0.5$'), join(images_dir, f'h_1e-10_ms{ms}_tau{tau}_s{s}.png'))

    #sw_1d(ms, tau, s)

    # Es = []
    # for tau in [0.01 / (2 ** i) for i in range(10)]:
    #     ts, m, E, me, Ee = sw_1d(ms, tau, s)
    #     Es.append(E[-1])
    # print(Es)

    # ts, m, E, me, Ee = sw_1d(ms, tau, s)

    # m_sc, E_sc, E1, E2 = check(ts)

    # plot_m(ts, me, m_sc, join(images_dir, f'ms{ms}_tau{tau}_s{s}_m.png'))
    # plot_E(ts, Ee, E_sc, join(images_dir, f'ms{ms}_tau{tau}_s{s}_E.png'))

    # plot_E1(ts, E1, join(images_dir, f'ms{ms}_tau{tau}_s{s}_E1.png'))
    # plot_E2(ts, E2, join(images_dir, f'ms{ms}_tau{tau}_s{s}_E2.png'))

    # ms = [100, 200, 400, 800, 1600]
    # Es = []
    # for m in ms:
    #     ts, m, E, me, Ee = sw_1d(m, tau, s)
    #     Es.append(Ee)
    # m_sc, E_sc, E1, E2 = check(ts)
    # Es.append(E_sc)
    # plot_Es(ts, Es, [*ms, 'analytic'], join(images_dir, f'Es.png'))
    # plot_Es_dif(ts, Es, [*ms, 'analytic'], join(images_dir, f'Es_dif.png'))



    # taus = [0.005, 0.0025, 0.00125, 0.000625]
    # ts, Es = [], []
    # E_dif = []
    # for tau in taus:
    #     t, m, E, me, Ee = sw_1d(ms, tau, s)
    #     ts.append(t)
    #     Es.append(Ee)
    #     #m_sc, E_sc, E1, E2 = check(t)
    #     E_dif.append(E - Ee)
    
    # t, m, E, me, Ee = sw_1d(ms, taus[-1], s)

    # ts.append(t)
    # Es.append(Ee)
    # plot_Es_tau(ts, Es, [*taus, 'analytic'], r'$E$', join(images_dir, f'Es_tau_1e-10.png'))
    # plot_Es_tau(ts[:-1], E_dif, taus, r'$E - E_{analytic}$', join(images_dir, f'Es_tau_dif_1e-10.png'))




    # print(np.allclose(me, m_sc), np.allclose(Ee, E_sc))
    # print(me)
    # print(m_sc)
    # print(Ee)
    # print(E_sc)
