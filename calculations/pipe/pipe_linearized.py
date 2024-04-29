from fenics import *
import numpy as np
from scipy.constants import pi, R


def calculate_pipe(mesh_size, tau):
    set_log_level(LogLevel.WARNING)

    P_in, P_out, m_in, m_out = [], [], [], []

    t_max = 48 * 3600
    L = 1e5
    T = 278

    D = 0.6
    A = pi * D**2 / 4

    eps = 0.000617
    Re = 6000
    f = (-2*np.log(eps/D/3.7 - 4.518/Re*np.log(6.9/Re + (eps/D/3.7)**1.11))) ** -2

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_left = 5e6
    P_0_right = 4.2e6

    rho = 0.73

    # t в часах
    m_out_expr = Expression('rho_out * (t < 2 ? 20*t+70 : (t < 10 ? 110 : (t < 12 ? -40*t+510 : (t < 22 ? 30 : 20*t-410))))', rho_out=rho, t=0, degree=1)

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

    Z = 1

    F = ((P-Pn)/(tau*Z*Rs*T) + m.dx(0)/A) * Pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/Pn).dx(0) + P.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*Pn)) * mt*dx
    a1, L1 = lhs(F), rhs(F)
    # print(a1)
    # print(L1)
    # exit()

    w = Function(W)

    def collect_data():
        P_in.append(w.sub(0)(0))
        P_out.append(w.sub(0)(L))
        m_in.append(w.sub(1)(0))
        m_out.append(w.sub(1)(L))

        print(f'Time {t:>7.5f}')

    t = 0
    w.assign(project(Expression(('(P_left + x[0]*(P_right - P_left)/L)', 'rho * 70'),
                                P_left=P_left, P_right=P_0_right, rho=rho, L=L, degree=1), W))
    collect_data()
    wn.assign(w)

    for t in ts[1:]:
        m_out_expr.t = t / 3600 if t / 3600 <= 24 else t / 3600 - 24
        solve(a1 == L1, w, bc)
        # solve(F == 0, w, bc, solver_parameters={"newton_solver": {
        #         'absolute_tolerance': 1e-5,
        #         'relative_tolerance': 1e-5,
        #         'maximum_iterations': 50,
        #         'relaxation_parameter': 1.0,
        #     }})
        collect_data()
        wn.assign(w)
    
    return map(np.array, (P_in, P_out, m_in, m_out))


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    P_in, P_out, m_in, m_out = calculate_pipe(mesh_size=100, tau=100)

    plt.figure()
    plt.plot(P_in)
    plt.title('P_in')

    plt.figure()
    plt.plot(P_out)
    plt.title('P_out')

    plt.figure()
    plt.plot(m_in)
    plt.title('m_in')

    plt.figure()
    plt.plot(m_out)
    plt.title('m_out')

    plt.show()
