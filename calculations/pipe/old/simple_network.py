from fenics import *
import numpy as np
from scipy.constants import pi, R


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau):
    t_max = 24 * 3600
    L_1 = 80000
    L_2 = 90000
    L_3 = 100000
    T = 278
    D = 0.6
    A = pi * D**2 / 4

    P_node_2, P_node_3 = [], []

    eps = 0.000617
    Re = 6000
    f = (-2*np.log(eps/D/3.7 - 4.518/Re*np.log(6.9/Re + (eps/D/3.7)**1.11))) ** -2
    f = 0.0105 # 0.003

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_left = 5e6
    q_right = 30
    rho = 0.7165
    Z = 1

    ts = np.arange(0, t_max+tau/2, tau)
    mesh1 = IntervalMesh(mesh_size, 0, L_1)
    mesh2 = IntervalMesh(mesh_size, 0, L_2)
    mesh3 = IntervalMesh(mesh_size, 0, L_3)

    P = FiniteElement('P', mesh1.ufl_cell(), 1)
    M = FiniteElement('P', mesh1.ufl_cell(), 1)
    W1 = FunctionSpace(mesh1, MixedElement([P, M]))
    W2 = FunctionSpace(mesh2, MixedElement([P, M]))
    W3 = FunctionSpace(mesh3, MixedElement([P, M]))
    
    bc1 = DirichletBC(W1.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')
    bc2 = DirichletBC(W2.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')

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

    w3n = Function(W3)
    pt, mt = TestFunctions(W3)
    p, m = TrialFunctions(W3)
    pn, mn = split(w3n)
    F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
    a3, L3 = lhs(F), rhs(F)
    w3 = Function(W3)

    def collect_data():
        if t >= 0:
            P_node_2.append(w3.sub(0)(0))
            P_node_3.append(w3.sub(0)(L_3))

    t = 0
    w1.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=rho*q_right, degree=1), W1))
    w2.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=rho*q_right, degree=1), W2))
    w3.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=rho*q_right, degree=1), W3))
    #collect_data()
    w1n.assign(w1)
    w2n.assign(w2)
    w3n.assign(w3)

    def m2():
        t_h = t / 3600
        q = None
        if t_h <= 4:
            q = 20 + 2.5*t_h
        elif t_h <= 12:
            q = 40 - 2.5*t_h
        elif t_h <= 20:
            q = 2.5*t_h - 20
        elif t_h <= 24:
            q = 80 - 2.5*t_h
        return rho * q

    def m3():
        t_h = t / 3600
        q = None
        if t_h <= 4:
            q = 40 + 2.5*t_h
        elif t_h <= 12:
            q = 60 - 2.5*t_h
        elif t_h <= 20:
            q = 2.5*t_h
        elif t_h <= 24:
            q = 100 - 2.5*t_h
        return rho * q

    
    # для установления течения, так как нет начальных условий
    print(tau, int(3600 // tau))
    for i in range(int(5 * 3600 // tau)):
        A1, b1 = assemble_system(a1, L1, bc1)
        A2, b2 = assemble_system(a2, L2, bc2)
        A3, b3 = assemble_system(a3, L3)

        A1, b1 = A1.array(), b1.get_local()
        A2, b2 = A2.array(), b2.get_local()
        A3, b3 = A3.array(), b3.get_local()

        n = b1.size
        n0 = (n + 1) * 3

        b0 = np.zeros(n0)
        b0[:n] = b3
        b0[n:2*n] = b2
        b0[2*n:3*n] = b1

        A0 = np.zeros((n0, n0))
        A0[:n, :n] = A3
        A0[n:2*n, n:2*n] = A2
        A0[2*n:3*n, 2*n:3*n] = A1

        v1, v2, v3 = n0-3, n0-2, n0-1

        # 1
        cur_n = 2 * n
        p1_end = cur_n
        m1_end = cur_n + 1
        A0[m1_end] = 0
        b0[m1_end] = 0
        A0[m1_end][m1_end] = 1
        A0[m1_end][v1] = -1

        # 2
        cur_n = n
        p2_end = cur_n
        m2_end = cur_n + 1
        A0[m2_end] = 0
        b0[m2_end] = 0
        A0[m2_end][m2_end] = 1
        A0[m2_end][v2] = -1

        # 3
        cur_n = 0
        p3_begin = cur_n + n - 2
        m3_begin = cur_n + n - 1
        p3_end = cur_n
        m3_end = cur_n + 1
        A0[p3_begin] = 0
        b0[p3_begin] = 0
        A0[p3_begin][p3_begin] = 1
        A0[p3_begin][v3] = -1

        A0[p3_end] = 0
        b0[p3_end] = 0
        A0[p3_end][p3_end] = 1
        A0[p3_end][p1_end] = -1

        #
        A0[v1][v1] = 1
        A0[v1][m3_end] = 1
        b0[v1] = m3()

        A0[v2][v2] = 1
        A0[v2][m3_begin] = -1
        b0[v2] = m2()

        A0[v3][v3] = 1
        A0[v3][p2_end] = -1

        res = np.linalg.solve(A0, b0)
        w3.vector()[:] = res[:n]
        w2.vector()[:] = res[n:2*n]
        w1.vector()[:] = res[2*n:3*n]

        w1n.assign(w1)
        w2n.assign(w2)
        w3n.assign(w3)

        #print(w1.sub(0)(L_1), w3.sub(0)(L_3), w2.sub(0)(L_2), w3.sub(0)(0))
    
    collect_data()

    for t in ts[1:]:
        A1, b1 = assemble_system(a1, L1, bc1)
        A2, b2 = assemble_system(a2, L2, bc2)
        A3, b3 = assemble_system(a3, L3)

        A1, b1 = A1.array(), b1.get_local()
        A2, b2 = A2.array(), b2.get_local()
        A3, b3 = A3.array(), b3.get_local()

        n = b1.size
        n0 = (n + 1) * 3

        b0 = np.zeros(n0)
        b0[:n] = b3
        b0[n:2*n] = b2
        b0[2*n:3*n] = b1

        A0 = np.zeros((n0, n0))
        A0[:n, :n] = A3
        A0[n:2*n, n:2*n] = A2
        A0[2*n:3*n, 2*n:3*n] = A1

        v1, v2, v3 = n0-3, n0-2, n0-1

        # 1
        cur_n = 2 * n
        p1_end = cur_n
        m1_end = cur_n + 1
        A0[m1_end] = 0
        b0[m1_end] = 0
        A0[m1_end][m1_end] = 1
        A0[m1_end][v1] = -1

        # 2
        cur_n = n
        p2_end = cur_n
        m2_end = cur_n + 1
        A0[m2_end] = 0
        b0[m2_end] = 0
        A0[m2_end][m2_end] = 1
        A0[m2_end][v2] = -1

        # 3
        cur_n = 0
        p3_begin = cur_n + n - 2
        m3_begin = cur_n + n - 1
        p3_end = cur_n
        m3_end = cur_n + 1
        A0[p3_begin] = 0
        b0[p3_begin] = 0
        A0[p3_begin][p3_begin] = 1
        A0[p3_begin][v3] = -1

        A0[p3_end] = 0
        b0[p3_end] = 0
        A0[p3_end][p3_end] = 1
        A0[p3_end][p1_end] = -1

        #
        A0[v1][v1] = 1
        A0[v1][m3_end] = 1
        b0[v1] = m3()

        A0[v2][v2] = 1
        A0[v2][m3_begin] = -1
        b0[v2] = m2()

        A0[v3][v3] = 1
        A0[v3][p2_end] = -1

        res = np.linalg.solve(A0, b0)
        w3.vector()[:] = res[:n]
        w2.vector()[:] = res[n:2*n]
        w1.vector()[:] = res[2*n:3*n]

        collect_data()

        w1n.assign(w1)
        w2n.assign(w2)
        w3n.assign(w3)

        #print(w1.sub(0)(L_1), w3.sub(0)(L_3), w2.sub(0)(L_2), w3.sub(0)(0))

    return w1.sub(0).compute_vertex_values(), w1.sub(1).compute_vertex_values(), \
        w2.sub(0).compute_vertex_values(), w2.sub(1).compute_vertex_values(), \
        w3.sub(0).compute_vertex_values(), w3.sub(1).compute_vertex_values(), \
        ts, P_node_2, P_node_3


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    p1, m1, p2, m2, p3, m3, ts, P_node_2, P_node_3 = calculate(mesh_size=50, tau=100) # 0.25 * 3600

    # fig, axs = plt.subplots(2, 3)
    # fig.set_size_inches(16, 9)
    # fig.set_dpi(300)

    # axs[0, 0].plot(p1)
    # axs[0, 0].set_title('P pipe 1')

    # axs[0, 1].plot(p2)
    # axs[0, 1].set_title('P pipe 2')

    # axs[0, 2].plot(p3)
    # axs[0, 2].set_title('P pipe 3')

    # axs[1, 0].plot(m1)
    # axs[1, 0].set_title('m pipe 1')

    # axs[1, 1].plot(m2)
    # axs[1, 1].set_title('m pipe 2')

    # axs[1, 2].plot(m3)
    # axs[1, 2].set_title('m pipe 3')

    fig, axs = plt.subplots(2)
    fig.set_size_inches(6, 6)
    fig.set_dpi(300)

    axs[0].plot(ts, P_node_2)
    axs[0].set_title('P_node_2')
    axs[0].grid()
    axs[0].set_xlim(ts[0], ts[-1])
    axs[0].set_ylim(4.7e6, 5e6)

    axs[1].plot(ts, P_node_3)
    axs[1].set_title('P_node_3')
    axs[1].grid()
    axs[1].set_xlim(ts[0], ts[-1])
    axs[1].set_ylim(4.7e6, 5e6)
    
    plt.show()
