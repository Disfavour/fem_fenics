from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, info=False):
    domain_size = 5

    u_= []
    m = []

    ts = np.arange(0, T+tau/2, tau)
    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)

    bc = DirichletBC(P, Constant(0), CompiledSubDomain('on_boundary && near(x[0], left)', left=-domain_size))

    u = Function(P)
    ut = TestFunction(P)
    un = Function(P)

    dt = (u - un)/tau
    us = sigma*u + (1-sigma)*un

    b = abs(u.dx(0)) ** gamma

    F = (dt + us.dx(0)) * ut*dx \
        + alfa*tau**2 * b * u.dx(0)*ut.dx(0) * dx(mesh) \
        + kappa * (u-un).dx(0) * ut*dx(mesh)
    # kappa*tau * us.dx(0) * ut*dx(mesh)

    def collect_data():
        m.append(assemble(u * dx))
        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f} m {m[-1]:>7.5f}')

    t = 0
    if ntest == 1:
        u.assign(project(Expression('x[0] >= -4.5 && x[0] <= -3.5 ? 1 : 0', degree=1), P))
    elif ntest == 2:
        u.assign(project(Expression('x[0] >= -4.5 && x[0] <= -3.5 ? sin((x[0] + 0.5) * pi) : 0', degree=degree), P))
    elif ntest == 3:
        u.assign(project(Expression('exp(-25*pow(x[0] + 4, 2))', degree=degree), P))
    collect_data()
    un.assign(u)

    for t in ts[1:]:
        solve(F == 0, u, bc, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-13,
            'relative_tolerance': 1e-13,
            'maximum_iterations': 50,
            'relaxation_parameter': 1.0,
        }})
        collect_data()
        un.assign(u)

    return map(np.array, (x, u_, ts, m))


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    strart_time = time.time()

    x, u, t, m = calculate(mesh_size=400, tau=0.01, degree=1, T=8.0, ts2store=[0.0, 4, 8.0], sigma=0.5, alfa=0, gamma=0, kappa=0.1, ntest=3, info=True)
    plt.figure()
    plt.plot(x, u.T)

    plt.figure()
    plt.plot(t, m)
    plt.xlim([0, 2])

    plt.show()

    # e^(-100800*(x-0.1)^(4))
    
    print(time.time() - strart_time)
    pass
