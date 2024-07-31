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

    # eps = 0.000617
    # Re = 6000
    # f = (-2*np.log(eps/D/3.7 - 4.518/Re*np.log(6.9/Re + (eps/D/3.7)**1.11))) ** -2
    f = 0.009

    S = 0.6
    M_air = 28.964917 / 1000
    M = S * M_air
    Rs = R / M

    P_in = 5e6
    rho = 0.73
    m_out = rho * 70

    ts = np.arange(0, t_max+tau/2, tau)
    mesh = IntervalMesh(mesh_size, 0, L)
    x = mesh.coordinates().flatten()

    P = FiniteElement('P', mesh.ufl_cell(), 1)
    M = FiniteElement('P', mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, MixedElement([P, M]))

    bc = [DirichletBC(W.sub(0), Constant(P_in), 'on_boundary && near(x[0], 0)'),
        DirichletBC(W.sub(1), Constant(m_out), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

    wn = Function(W)
    pt, mt = TestFunctions(W)
    p, m = TrialFunctions(W)
    pn, mn = split(wn)
    
    F = ((p-pn)/(tau*Rs*T) + m.dx(0)/A) * pt*dx \
        + ((m-mn)/(tau*A) + Rs*T/A**2*(m*mn/pn).dx(0) + p.dx(0) + Rs*T*f*m*abs(mn)/(2*D*A**2*pn)) * mt*dx
    a1, L1 = lhs(F), rhs(F)

    w = Function(W)

    t = 0
    w.assign(project(Expression(('p', 'm'), p=P_in, m=m_out, degree=1), W))
    wn.assign(w)

    for t in ts[1:]:
        solve(a1 == L1, w, bc)
        wn.assign(w)

    return x, w.sub(0).compute_vertex_values(), stationary_analytic.get_exact(x, D, A, Rs, T, f, m_out, P_in)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, p_numerical, p_analytic = calculate(mesh_size=100, tau=4 * 3600, t_max=4 * 3600)

    print(np.allclose(p_numerical, p_analytic))

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(x, p_analytic, 'b')
    plt.plot(x, p_numerical, '--r')
    plt.legend(['analytic', 'numerical'])
    plt.title('P')

    plt.show()
