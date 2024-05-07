from fenics import *
import numpy as np
from scipy.constants import pi, R


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau):
    t_max = 12 * 3600
    L = 36300
    T = 278
    D = 1.422
    A = pi * D**2 / 4

    P_nodes = [[] for i in range(6)]
    m_in = [[] for i in range(7)]
    m_mid = [[] for i in range(7)]
    m_out = [[] for i in range(7)]

    eps = 0.000015
    Re = 5000
    f = (-2*np.log(eps/D/3.7 - 4.518/Re*np.log(6.9/Re + (eps/D/3.7)**1.11))) ** -2

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_in = 8e6
    Z = 1

    # t в часах
    m_out_expr = Expression('t < 1 ? 100 : (t < 3 ? 50*t + 50 : (t < 6 ? 200 : (t < 8 ? 560 - 60*t : 80)))', t=0, degree=1)

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)

    P = FiniteElement('P', mesh.ufl_cell(), 1)
    M = FiniteElement('P', mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, MixedElement([P, M]))

    bc0 = DirichletBC(W.sub(0), Constant(P_in), 'on_boundary && near(x[0], 0)')
    bc1, bc2, bc3, bc4, bc5 = [], [], [], [], []
    bc6 = DirichletBC(W.sub(1), m_out_expr, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))
    bc_s = [bc0, bc1, bc2, bc3, bc4, bc5, bc6]
    npipes = len(bc_s)

    w_s, wn_s = [], []
    a_s, L_s = [], []
    for i in range(npipes):
        p, m = TrialFunctions(W)
        pt, mt = TestFunctions(W)
        w_s.append(Function(W))
        wn_s.append(Function(W))
        pn, mn = split(wn_s[-1])
        F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
            + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
        a_s.append(lhs(F))
        L_s.append(rhs(F))

    def collect_data():
        print(f'Time {t:>7.5f}')
        P_nodes[0].append(w_s[0].sub(0)(0))
        P_nodes[1].append(w_s[1].sub(0)(0))
        P_nodes[2].append(w_s[4].sub(0)(0))
        P_nodes[3].append(w_s[5].sub(0)(0))
        P_nodes[4].append(w_s[6].sub(0)(0))
        P_nodes[5].append(w_s[6].sub(0)(L))

        for i in range(npipes):
            m_in[i].append(w_s[i].sub(1)(0))
            m_mid[i].append(w_s[i].sub(1)(L/2))
            m_out[i].append(w_s[i].sub(1)(L))
    
    t = 0
    for w, wn in zip(w_s, wn_s):
        w.assign(project(Expression(('P', 'm'), P=P_in, m=100, degree=1), W))
        wn.assign(w)
    
    def iteration():
        A_s, b_s = [], []
        for a, L, bc in zip(a_s, L_s, bc_s):
            A, b = assemble_system(a, L, bc)
            A_s.append(A.array())
            b_s.append(b.get_local())

        n = b_s[0].size
        N = npipes * (n + 2) - 2    # 12

        b0 = np.zeros(N)
        for i, b in enumerate(reversed(b_s)):
            b0[i*n:(i+1)*n] = b
        
        A0 = np.zeros((N, N))
        for i, A in enumerate(reversed(A_s)):
            A0[i*n:(i+1)*n, i*n:(i+1)*n] = A
        
        v = []
        for i in range(2*npipes - 2):
            v.append(N - 1 - i)
        v.reverse()

        # 6
        m = 0
        p6_in = m + n - 2
        A0[p6_in] = 0
        b0[p6_in] = 0
        A0[p6_in][p6_in] = 1
        A0[p6_in][v[0]] = -1
        p6_in = v[0]
        m6_in = m + n - 1

        # 5
        m += n

        p5_in = m + n - 2
        A0[p5_in] = 0
        b0[p5_in] = 0
        A0[p5_in][p5_in] = 1
        A0[p5_in][v[1]] = -1
        p5_in = v[1]
        m5_in = m + n - 1

        m5_out = m + 1
        A0[m5_out] = 0
        b0[m5_out] = 0
        A0[m5_out][m5_out] = 1
        A0[m5_out][v[2]] = -1
        m5_out = v[2]
        p5_out = m

        # 4
        m += n

        p4_in = m + n - 2
        A0[p4_in] = 0
        b0[p4_in] = 0
        A0[p4_in][p4_in] = 1
        A0[p4_in][v[3]] = -1
        p4_in = v[3]
        m4_in = m + n - 1

        m4_out = m + 1
        A0[m4_out] = 0
        b0[m4_out] = 0
        A0[m4_out][m4_out] = 1
        A0[m4_out][v[4]] = -1
        m4_out = v[4]
        p4_out = m

        # 3
        m += n

        p3_in = m + n - 2
        A0[p3_in] = 0
        b0[p3_in] = 0
        A0[p3_in][p3_in] = 1
        A0[p3_in][v[5]] = -1
        p3_in = v[5]
        m3_in = m + n - 1

        m3_out = m + 1
        A0[m3_out] = 0
        b0[m3_out] = 0
        A0[m3_out][m3_out] = 1
        A0[m3_out][v[6]] = -1
        m3_out = v[6]
        p3_out = m

        # 2
        m += n

        p2_in = m + n - 2
        A0[p2_in] = 0
        b0[p2_in] = 0
        A0[p2_in][p2_in] = 1
        A0[p2_in][v[7]] = -1
        p2_in = v[7]
        m2_in = m + n - 1

        m2_out = m + 1
        A0[m2_out] = 0
        b0[m2_out] = 0
        A0[m2_out][m2_out] = 1
        A0[m2_out][v[8]] = -1
        m2_out = v[8]
        p2_out = m

        # 1
        m += n

        p1_in = m + n - 2
        A0[p1_in] = 0
        b0[p1_in] = 0
        A0[p1_in][p1_in] = 1
        A0[p1_in][v[9]] = -1
        p1_in = v[9]
        m1_in = m + n - 1

        m1_out = m + 1
        A0[m1_out] = 0
        b0[m1_out] = 0
        A0[m1_out][m1_out] = 1
        A0[m1_out][v[10]] = -1
        m1_out = v[10]
        p1_out = m

        # 0
        m += n
        m0_out = m + 1
        A0[m0_out] = 0
        b0[m0_out] = 0
        A0[m0_out][m0_out] = 1
        A0[m0_out][v[11]] = -1
        m0_out = v[11]
        p0_out = m

        # node 1
        A0[v[0]][p0_out] = 1
        A0[v[0]][p1_in] = -1

        A0[v[1]][p0_out] = 1
        A0[v[1]][p2_in] = -1

        A0[v[2]][m0_out] = 1
        A0[v[2]][m1_in] = -1
        A0[v[2]][m2_in] = -1

        # node 2
        A0[v[3]][p1_out] = 1
        A0[v[3]][p3_in] = -1

        A0[v[4]][p1_out] = 1
        A0[v[4]][p4_in] = -1

        A0[v[5]][m1_out] = 1
        A0[v[5]][m3_in] = -1
        A0[v[5]][m4_in] = -1

        # node 3
        A0[v[6]][p2_out] = 1
        A0[v[6]][p3_out] = -1

        A0[v[7]][p2_out] = 1
        A0[v[7]][p5_in] = -1

        A0[v[8]][m2_out] = 1
        A0[v[8]][m3_out] = 1
        A0[v[8]][m5_in] = -1

        # node 4
        A0[v[9]][p4_out] = 1
        A0[v[9]][p5_out] = -1

        A0[v[10]][p4_out] = 1
        A0[v[10]][p6_in] = -1

        A0[v[11]][m4_out] = 1
        A0[v[11]][m5_out] = 1
        A0[v[11]][m6_in] = -1

        res = np.linalg.solve(A0, b0)
        for i, w in enumerate(reversed(w_s)):
            w.vector()[:] = res[i*n:(i+1)*n]
        
        for w, wn in zip(w_s, wn_s):
            wn.assign(w)
    
    # установившееся течение - это начальные условия
    for i in range(int(3600 // tau)):
        iteration()
    
    collect_data()

    for t in ts[1:]:
        m_out_expr.t = t / 3600
        iteration()
        collect_data()

    return ts, P_nodes, m_in, m_mid, m_out


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    t, P_nodes, m_in, m_mid, m_out = calculate(mesh_size=50, tau=100) # 0.25 * 3600

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    for i, p in enumerate(P_nodes):
        plt.plot(t, p, label=i)
    plt.legend()
    plt.grid()

    for ms in [m_out]:
        plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
        for i, m in enumerate(ms):
            plt.plot(t, m, label=i)
        plt.legend()
        plt.grid()
    
    
    plt.show()
