from fenics import *
import numpy as np
from scipy.constants import pi, R
from bc import BC_periodic


def calculate_pipe(mesh_size, tau):
    set_log_level(LogLevel.WARNING)

    P_in, P_out, m_in, m_out = [], [], [], []

    t_max = 48 * 3600
    L = 1e5
    T = 278

    D = 0.6
    A = pi * D**2 / 4

    rho = 0.73

    f = 0.009 #0.00875

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_left = 5e6

    m_out_expr = BC_periodic(rho)

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)

    P_ = FiniteElement('P', mesh.ufl_cell(), 1)
    m_ = FiniteElement('P', mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, MixedElement([P_, m_]))

    bc = [DirichletBC(W.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)'),
        DirichletBC(W.sub(1), m_out_expr, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

    wn = Function(W)
    Pt, mt = TestFunctions(W)
    P, m = TrialFunctions(W)
    Pn, mn = split(wn)

    F = ((P-Pn)/(tau*Rs*T) + m.dx(0)/A) * Pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(m*mn/Pn).dx(0) + P.dx(0) + Rs*T*f*m*abs(mn)/(2*D*A**2*Pn)) * mt*dx
    a1, L1 = lhs(F), rhs(F)

    w = Function(W)

    def collect_data():
        print(f'Time {t:>7.5f}')
        P_in.append(w.sub(0)(0))
        P_out.append(w.sub(0)(L))
        m_in.append(w.sub(1)(0))
        m_out.append(w.sub(1)(L))

    t = 0
    w.assign(project(Expression(('P', 'm'), P=P_left, m=rho*70, degree=0), W))

    # установившееся течение - это начальные условия
    while not np.allclose(w.vector().get_local(), wn.vector().get_local()):
        wn.assign(w)
        solve(a1 == L1, w, bc)
        t += tau
        print(f'Time ustanovlenie {t:>7.5f}')
    
    t = 0
    collect_data()
    wn.assign(w)

    for t in ts[1:]:
        m_out_expr.update_t(t)
        solve(a1 == L1, w, bc)
        collect_data()
        wn.assign(w)
    
    return map(np.array, (ts, P_in, P_out, m_in, m_out))


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    t, P_in, P_out, m_in, m_out = calculate_pipe(mesh_size=100, tau=100)

    # plt.figure()
    # plt.plot(t, P_in)
    # plt.title('P_in')
    # plt.grid()

    plt.figure()
    plt.plot(t, P_out)
    plt.title('P_out')
    plt.grid()
    plt.ylim(2.5e6, 5e6)
    plt.xlim(t[0], t[-1])

    # plt.figure()
    # plt.plot(t, m_in)
    # plt.title('m_in')
    # plt.grid()

    # plt.figure()
    # plt.plot(t, m_out)
    # plt.title('m_out')
    # plt.grid()

    plt.show()
