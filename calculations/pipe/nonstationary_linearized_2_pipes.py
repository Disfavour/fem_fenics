from fenics import *
import numpy as np
from scipy.constants import pi, R
import stationary_analytic

set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, t_max):
    L = 1e5 // 2
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
    rho = 0.73
    m_right = rho * 70
    Z = 1

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)
    x = mesh.coordinates().flatten()

    P = FiniteElement('P', mesh.ufl_cell(), 1)
    M = FiniteElement('P', mesh.ufl_cell(), 1)
    W1 = FunctionSpace(mesh, MixedElement([P, M]))
    W2 = FunctionSpace(mesh, MixedElement([P, M]))
    
    bc1 = DirichletBC(W1.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')
    bc2 = DirichletBC(W2.sub(1), Constant(m_right), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))

    w1n = Function(W1)
    pt, mt = TestFunctions(W1)
    p, m = TrialFunctions(W1)
    pn, mn = split(w1n)
    F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
    a1, L1 = lhs(F), rhs(F)
    w1 = Function(W1)

    w2n = Function(W2)
    pt, mt = TestFunctions(W2)
    p, m = TrialFunctions(W2)
    pn, mn = split(w2n)
    F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
    a2, L2 = lhs(F), rhs(F)
    w2 = Function(W2)

    t = 0
    w1.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W1))
    w2.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W2))
    w1n.assign(w1)
    w2n.assign(w2)

    for t in ts[1:]:
        A1, b1 = assemble_system(a1, L1, bc1)
        A2, b2 = assemble_system(a2, L2, bc2)

        A1, b1 = A1.array(), b1.get_local()
        A2, b2 = A2.array(), b2.get_local()

        n = b1.size
        n0 = (n + 1) * 2

        b0 = np.zeros(n0)
        b0[:n] = b2
        b0[n:2*n] = b1

        A0 = np.zeros((n0, n0))
        A0[:n, :n] = A2
        A0[n:2*n, n:2*n] = A1

        # P на начале 2-й трубы
        p2_left = np.zeros(n0)
        p2_left[n-2] = 1
        p2_left[n0-2] = -1
        A0[n-2] = p2_left
        b0[n-2] = 0
        # равно P на конце 1-й трубы
        p2_left = np.zeros(n0)
        p2_left[n0-2] = 1
        p2_left[n] = -1
        A0[n0-2] = p2_left
        b0[n0-2] = 0

        # m на конце 1-й трубы
        m1_right = np.zeros(n0)
        m1_right[n+1] = 1
        m1_right[n0-1] = -1
        A0[n+1] = m1_right
        b0[n+1] = 0
        # равно m в начале 2-й трубы
        m1_right = np.zeros(n0)
        m1_right[n0-1] = 1
        m1_right[n-1] = -1
        A0[n0-1] = m1_right
        b0[n0-1] = 0

        res = np.linalg.solve(A0, b0)
        w2.vector()[:] = res[:n]
        w1.vector()[:] = res[n:2*n]

        w2n.assign(w2)
        w1n.assign(w1)

    return x, w1.sub(0).compute_vertex_values(), w1.sub(1).compute_vertex_values(), w2.sub(0).compute_vertex_values(), w2.sub(1).compute_vertex_values()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, p1, m1, p2, m2 = calculate(mesh_size=10, tau=0.25 * 3600, t_max=8 * 3600)

    #plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    fig, axs = plt.subplots(2, 2)
    fig.set_size_inches(16, 9)
    fig.set_dpi(300)

    axs[0, 0].plot(x, p1)
    axs[0, 0].set_title('P pipe 1')

    axs[0, 1].plot(x, p2)
    axs[0, 1].set_title('P pipe 2')

    axs[1, 0].plot(x, m1)
    axs[1, 0].set_title('m pipe 1')

    axs[1, 1].plot(x, m2)
    axs[1, 1].set_title('m pipe 2')
    
    plt.show()

    # plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    # plt.plot(x, p1)
    # plt.title("P pipe 1")

    # plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    # plt.plot(x, m1)
    # plt.title("m pipe 1")

    # plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    # plt.plot(x, p2)
    # plt.title("P pipe 2")

    # plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    # plt.plot(x, m2)
    # plt.title("m pipe 2")

    # plt.show()
