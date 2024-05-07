from fenics import *
import numpy as np
from scipy.constants import pi, R
import stationary_analytic

set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, t_max):
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
    rho = 0.73
    m_right = rho * 70
    Z = 1

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)
    x = mesh.coordinates().flatten()
    mesh1 = IntervalMesh((mesh_size - 1) // 2, 0, L/2)
    mesh2 = IntervalMesh((mesh_size - 1) // 2, L/2, L)

    P = FiniteElement('P', mesh.ufl_cell(), 1)
    M = FiniteElement('P', mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, MixedElement([P, M]))
    W1 = FunctionSpace(mesh1, MixedElement([P, M]))
    W2 = FunctionSpace(mesh2, MixedElement([P, M]))

    bc = [DirichletBC(W.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)'),
        DirichletBC(W.sub(1), Constant(m_right), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]
    bc1 = DirichletBC(W1.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)')
    bc2 = DirichletBC(W2.sub(1), Constant(m_right), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))

    wn = Function(W)
    pt, mt = TestFunctions(W)
    p, m = TrialFunctions(W)
    pn, mn = split(wn)
    F = ((p-pn)/(tau*Z*Rs*T) + m.dx(0)/A) * pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(Z*m*mn/pn).dx(0) + p.dx(0) + Z*Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
    a0, L0 = lhs(F), rhs(F)
    w = Function(W)

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
    #w.assign(project(Expression(('P_left-30*x[0]', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W))
    w.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W))
    w1.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W1))
    w2.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W2))

    #w.vector()[dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()]] = exact
    wn.assign(w)
    w1n.assign(w1)
    w2n.assign(w2)

    # A = assemble(a1)
    # b = assemble(L1)
    # bc.apply(A, b)

    A, b = assemble_system(a0, L0, bc)
    print(A.array(), b.get_local(), sep='\n')
    A1, b1 = assemble_system(a1, L1, bc1)
    print(A1.array(), b1.get_local(), sep='\n')
    A2, b2 = assemble_system(a2, L2, bc2)
    print(A2.array(), b2.get_local(), sep='\n')

    res1 = np.linalg.solve(A1.array(), b1.get_local())
    res2 = np.linalg.solve(A2.array(), b2.get_local())


    #print(as_backend_type(A).mat().getValue(1, 1))
    #solve(A.array(), w.vector(), b.get_local())

    for t in ts[1:]:
        #solve(a1 == L1, w, bc)
        A1, b1 = assemble_system(a1, L1, bc1)
        A2, b2 = assemble_system(a2, L2, bc2)

        A0 = np.zeros_like(A.array())
        n = A0.shape[0]
        m = n // 2

        A0[:m, :m] = A2.array()
        A0[m:, m:] = A1.array()

        b0 = np.zeros_like(b.get_local())
        b0[:m] = b2.get_local()
        b0[m:] = b1.get_local()

        # link
        #b0[m-2:m] = 0

        link_P = np.zeros_like(A0[0])
        link_P[m-2] = 1
        link_P[m] = -1
        A0[m-2] = link_P
        b0[m-2] = 0
        # вторая труба получает граничное условие P слева (в начале)

        link_m = np.zeros_like(A0[0])
        link_m[m-1] = 1
        link_m[m+1] = -1
        A0[m+1] = link_m
        b0[m+1] = 0
        # первая труба получает граничное условие m справа (на конце)

        # print(A0)
        # print(b0)
        # exit()

        res = np.linalg.solve(A0, b0)
        w2.vector()[:] = res[:m]
        w1.vector()[:] = res[m:]

        w2n.assign(w2)
        w1n.assign(w1)

        #wn.assign(w)
        #print(w.vector()[:])
        #exit()
    w.vector()[:] = np.concatenate((w2n.vector()[:], w1n.vector()[:]), axis=None)
    #w.vector()[:m] = w2n.vector()[:]
    #w.vector()[m:] = w1n.vector()[:]
    import matplotlib.pyplot as plt

    plt.figure()
    plot(w.sub(0))
    plt.figure()
    plot(w.sub(1))
    plt.show()

    exit()

    return x, w.sub(0).compute_vertex_values(), stationary_analytic.get_exact(x, D, A, Rs, T, f, m_right, P_left), w.sub(1).compute_vertex_values(), np.full_like(x, m_right)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, p_numerical, p_analytic, m_numerical, m_analytic = calculate(mesh_size=11, tau=0.25 * 3600, t_max=8 * 3600)

    print(np.allclose(p_numerical, p_analytic))

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(x, p_analytic, 'b')
    plt.plot(x, p_numerical, '--r')
    plt.legend(['analytic', 'numerical'])
    plt.title('P')

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(x, m_numerical, 'b')
    plt.plot(x, m_analytic, '--r')
    plt.legend(['analytic', 'numerical'])
    plt.title('m')

    plt.show()
