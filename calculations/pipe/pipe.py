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

    f = 0.00875

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_left = 5e6

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

    w, wn = Function(W), Function(W)
    Pt, mt = TestFunctions(W)
    P, m = split(w)
    Pn, mn = split(wn)

    Z = 1

    F = ((P-Pn)/(tau*Z*Rs*T) + m.dx(0)/A) * Pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m**2/P).dx(0) + P.dx(0) + Z*Rs*T*f*m*abs(m)/(2*D*A**2*P)) * mt*dx

    def collect_data():
        print(f'Time {t:>7.5f}')
        P_in.append(w.sub(0)(0))
        P_out.append(w.sub(0)(L))
        m_in.append(w.sub(1)(0))
        m_out.append(w.sub(1)(L))

    t = 0
    w.assign(project(Expression(('P', 'm'), P=P_left, m=rho*70, degree=0), W))
    wn.assign(w)

    # установившееся течение - это начальные условия
    for i in range(int(5*3600 // tau)):
        solve(F == 0, w, bc, solver_parameters={"newton_solver": {
                'absolute_tolerance': 1e-5,
                'relative_tolerance': 1e-5,
                'maximum_iterations': 50,
                'relaxation_parameter': 1.0,
            }})
        wn.assign(w)
    
    collect_data()

    for t in ts[1:]:
        m_out_expr.t = t / 3600 if t / 3600 <= 24 else t / 3600 - 24
        solve(F == 0, w, bc, solver_parameters={"newton_solver": {
                'absolute_tolerance': 1e-5,
                'relative_tolerance': 1e-5,
                'maximum_iterations': 50,
                'relaxation_parameter': 1.0,
            }})
        collect_data()
        wn.assign(w)
    
    return map(np.array, (ts, P_in, P_out, m_in, m_out))


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    t, P_in, P_out, m_in, m_out = calculate_pipe(mesh_size=100, tau=100)

    plt.figure()
    plt.plot(t, P_out)
    plt.title('P_out')
    plt.grid()
    plt.ylim(2.5e6, 5e6)
    plt.xlim(t[0], t[-1])

    plt.show()
