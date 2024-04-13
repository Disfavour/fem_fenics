from fenics import *
import numpy as np
from scipy.constants import pi, R
import stationary_analytic


set_log_level(LogLevel.WARNING)


def calculate(mesh_size):
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

    m = rho * 70

    mesh = IntervalMesh(mesh_size, 0, L)
    x = mesh.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', 1)

    bc = DirichletBC(P, Constant(P_left), 'on_boundary && near(x[0], 0)')

    p = Function(P)
    pt = TestFunction(P)

    F = (p.dx(0) + Rs*T*m**2/A**2 * (p**-1).dx(0) + Rs*T*f*m**2/(2*D*A**2*p)) * pt*dx

    p.assign(project(Constant(P_left), P))

    solve(F == 0, p, bc)
    mesh.coordinates

    return x, p.compute_vertex_values(), stationary_analytic.get_exact(x, D, A, Rs, T, f, m, P_left)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, p_numerical, p_analytic = calculate(mesh_size=100)

    print(np.allclose(p_numerical, p_analytic))

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(x, p_analytic, 'b')
    plt.plot(x, p_numerical, '--r')
    plt.legend(['analytic', 'numerical'])
    plt.show()
