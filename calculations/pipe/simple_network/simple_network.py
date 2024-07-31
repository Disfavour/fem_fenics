from fenics import *
import numpy as np
from scipy.constants import pi, R
from bc import get_bcs


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

    f = 0.009    # 0.0105

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_left = 5e6
    q_right = 40
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
    Ws = [W1, W2, W3]
    
    bc0 = DirichletBC(W1.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')
    bc1 = DirichletBC(W2.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')
    bc2 = []
    bc = [bc0, bc1, bc2]
    npipes = len(bc)

    w, wn = [], []
    a, L = [], []
    for W in Ws:
        p, m = TrialFunctions(W)
        pt, mt = TestFunctions(W)
        w.append(Function(W))
        wn.append(Function(W))
        pn, mn = split(wn[-1])
        F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
            + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
        a.append(lhs(F))
        L.append(rhs(F))
    
    def collect_data():
        print(f'Time {t:>7.5f}')
        P_node_2.append(w[2].sub(0)(0))
        P_node_3.append(w[2].sub(0)(L_3))
    
    q2, q3 = get_bcs()

    def calculate_w():
        A, b = [], []
        for a1, L1, bc1 in zip(a, L, bc):
            A1, b1 = assemble_system(a1, L1, bc1)
            A.append(A1.array())
            b.append(b1.get_local())

        n = b[0].size
        N = npipes*n + 4

        b0 = np.zeros(N)
        for i, b1 in enumerate(reversed(b)):
            b0[i*n:(i+1)*n] = b1
        
        A0 = np.zeros((N, N))
        for i, A1 in enumerate(reversed(A)):
            A0[i*n:(i+1)*n, i*n:(i+1)*n] = A1
        
        v1, v2, v3, v4 = N-4, N-3, N-2, N-1

        # 1
        cur_n = 2 * n
        m1_end = cur_n + 1
        A0[m1_end] = 0
        b0[m1_end] = 0
        A0[m1_end][m1_end] = 1
        A0[m1_end][v1] = -1

        p1_end = cur_n
        m1_end = v1

        # 2
        cur_n = n
        m2_end = cur_n + 1
        A0[m2_end] = 0
        b0[m2_end] = 0
        A0[m2_end][m2_end] = 1
        A0[m2_end][v2] = -1

        p2_end = cur_n
        m2_end = v2

        # 3
        cur_n = 0

        p3_begin = cur_n + n - 2
        A0[p3_begin] = 0
        b0[p3_begin] = 0
        A0[p3_begin][p3_begin] = 1
        A0[p3_begin][v3] = -1

        p3_begin = v3
        m3_begin = cur_n + n - 1

        m3_end = cur_n + 1
        A0[m3_end] = 0
        b0[m3_end] = 0
        A0[m3_end][m3_end] = 1
        A0[m3_end][v4] = -1

        p3_end = cur_n
        m3_end = v4

        # node 2
        A0[v1][p2_end] = 1
        A0[v1][p3_begin] = -1

        A0[v2][m2_end] = 1
        A0[v2][m3_begin] = -1
        b0[v2] = rho * q2(t)

        # node 3
        A0[v3][p1_end] = 1
        A0[v3][p3_end] = -1

        A0[v4][m1_end] = 1
        A0[v4][m3_end] = 1
        b0[v4] = rho * q3(t)

        res = np.linalg.solve(A0, b0)
        for i, w1 in enumerate(reversed(w)):
            w1.vector()[:] = res[i*n:(i+1)*n]
    
    t = 0
    t_ust = 0
    for w_cur, wn_cur, W in zip(w, wn, Ws):
        w_cur.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=rho*q_right, degree=1), W))
    
    # установившееся течение - это начальные условия
    while not np.allclose(np.array([w_cur.vector().get_local() for w_cur in w]), np.array([wn_cur.vector().get_local() for wn_cur in wn])):
        for w_cur, wn_cur in zip(w, wn):
            wn_cur.assign(w_cur)
        
        calculate_w()

        t_ust += tau
        print(f'Time ustanovlenie {t_ust:>7.5f}')
    
    collect_data()
    for w_cur, wn_cur in zip(w, wn):
            wn_cur.assign(w_cur)

    for t in ts[1:]:
        calculate_w()
        
        collect_data()

        for w_cur, wn_cur in zip(w, wn):
            wn_cur.assign(w_cur)

    return ts, P_node_2, P_node_3


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from os.path import join

    taus = [i * 3600 / 4 for i in (1, 2, 4)]
    print(taus)

    lines = ['-', '--', ':', '-.']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    for tau, l in zip(taus, lines):
        ts, P_node_2, P_node_3 = calculate(mesh_size=100, tau=tau)
        plt.plot(ts, P_node_2, colors[0]+l)
        plt.plot(ts, P_node_3, colors[1]+l)

    plt.xticks(range(0, 25*3600, 4*3600))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$p$')
    plt.xlim(ts[0], ts[-1])
    plt.legend(('2', '3'))
    plt.grid()
    plt.savefig(join('images', 'pipes', 'simple.pdf'), transparent=True)
