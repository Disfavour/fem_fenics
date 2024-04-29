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
    m_right2 = rho * 40
    m_right3 = rho * 30
    Z = 1

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)
    x = mesh.coordinates().flatten()

    P = FiniteElement('P', mesh.ufl_cell(), 1)
    M = FiniteElement('P', mesh.ufl_cell(), 1)
    W1 = FunctionSpace(mesh, MixedElement([P, M]))
    W2 = FunctionSpace(mesh, MixedElement([P, M]))
    W3 = FunctionSpace(mesh, MixedElement([P, M]))

    m1 = Expression('m', m=m_right, degree=0)
    p23 = Expression('p', p=P_left, degree=0)
    bc1 = [
        DirichletBC(W1.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)'),
        DirichletBC(W1.sub(1), m1, CompiledSubDomain('on_boundary && near(x[0], L)', L=L))
    ]
    bc2 = [
        DirichletBC(W2.sub(0), p23, 'on_boundary && near(x[0], 0)'),
        DirichletBC(W2.sub(1), Constant(m_right2), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))
    ]
    bc3 = [
        DirichletBC(W3.sub(0), p23, 'on_boundary && near(x[0], 0)'),
        DirichletBC(W3.sub(1), Constant(m_right3), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))
    ]

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

    t = 0
    # w1.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W1))
    # w2.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W2))
    # w3.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W3))

    w1.assign(project(Expression(('P_left-6.8*x[0]', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W1))
    w2.assign(project(Expression(('-2.4*x[0] + 4660000', 'm_right'), P_left=P_left, m_right=m_right2, degree=1), W2))
    w3.assign(project(Expression(('-1.4*x[0] + 4660000', 'm_right'), P_left=P_left, m_right=m_right3, degree=1), W3))

    w1n.assign(w1)
    w2n.assign(w2)
    w3n.assign(w3)

    for t in ts[1:]:
        m1.m = w2.sub(1)(0) + w3.sub(1)(0)
        p23.p = w1.sub(0)(L)

        solve(a1 == L1, w1, bc1)

        #p23.p = w1.sub(0)(L)
        #print(t, m1.m, p23.p, w2.sub(1)(0), w3.sub(1)(0))

        solve(a2 == L2, w2, bc2)
        solve(a3 == L3, w3, bc3)

        w1n.assign(w1)
        w2n.assign(w2)
        w3n.assign(w3)

    return x, w1.sub(0).compute_vertex_values(), w1.sub(1).compute_vertex_values(), \
        w2.sub(0).compute_vertex_values(), w2.sub(1).compute_vertex_values(), \
        w3.sub(0).compute_vertex_values(), w3.sub(1).compute_vertex_values()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, p1, m1, p2, m2, p3, m3 = calculate(mesh_size=10, tau=0.25 * 3600, t_max=8 * 3600)

    fig, axs = plt.subplots(2, 3)
    fig.set_size_inches(16, 9)
    fig.set_dpi(300)

    axs[0, 0].plot(x, p1)
    axs[0, 0].set_title('P pipe 1')

    axs[0, 1].plot(x, p2)
    axs[0, 1].set_title('P pipe 2')

    axs[0, 2].plot(x, p3)
    axs[0, 2].set_title('P pipe 3')

    axs[1, 0].plot(x, m1)
    axs[1, 0].set_title('m pipe 1')

    axs[1, 1].plot(x, m2)
    axs[1, 1].set_title('m pipe 2')

    axs[1, 2].plot(x, m3)
    axs[1, 2].set_title('m pipe 3')
    
    plt.show()
