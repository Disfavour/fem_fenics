from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, info=False):
    if ntest == 1:
        a, b, x0 = 0.5, 1.5, 0.2
    elif ntest == 2:
        a, b, x0 = 1.5, 0.5, 0.2

    u_, ue_, err = [], [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = UnitIntervalMesh(mesh_size)
    x = mesh.coordinates().flatten()
    mesh_exact = UnitIntervalMesh(10000)
    xe = mesh_exact.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)
    PE = FunctionSpace(mesh_exact, 'P', 1)

    bc = DirichletBC(P, Constant(a), CompiledSubDomain('on_boundary && near(x[0], 0)'))

    u, un = Function(P), Function(P)
    ut = TestFunction(P)
    ue = Function(PE)

    dt = (u - un)/tau
    us = sigma*u + (1-sigma)*un

    b_ = abs(us.dx(0)) ** gamma

    F = (dt + us*us.dx(0)) * ut*dx \
        + alfa*tau**2 * b_ * us.dx(0)*ut.dx(0) * dx(mesh) \
        + kappa * (u-un).dx(0) * ut*dx(mesh)
    
    def exact_rarefaction(t):
        xs = xe.copy()
        if t > 0:
            return (xs <= x0 + a*t)*a + (xs > x0 + a*t)*(xs <= x0 + b*t)*(xs - x0)/t + (xs > x0 + b*t) * b
        else:
            return (xs <= 0.2)*a + (xs > 0.2)*b
    
    def exact_shock(t):
        xs = xe.copy()
        D = (a + b) / 2
        if t > 0:
            return (xs <= x0 + D*t)*a + (xs > x0 + D*t) * b
        else:
            return (xs <= 0.2) * a + (xs > 0.2) * b
    
    if ntest == 1:
        exact = exact_rarefaction
    elif ntest == 2:
        exact = exact_shock

    def collect_data():
        ue.vector()[dof_to_vertex_map(PE)[PE.dofmap().dofs()]] = exact(t)

        err.append(errornorm(u, ue, norm_type='L2', degree_rise=0))

        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())
            ue_.append(ue.compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f}')

    t = 0
    u.assign(project(Expression('x[0] <= x0 ? a : b', x0=x0, a=a, b=b, degree=1), P))
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

    return map(np.array, (x, u_, xe, ue_, ts, err))


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    strart_time = time.time()

    x, u_, xe, ue_, ts, err = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.15, 0.3], sigma=0.5, alfa=1000, gamma=0, kappa=0.0, ntest=2, info=True)
    plt.plot(x, u_.T)
    plt.plot(xe, ue_.T, '--')

    plt.figure()
    plt.plot(ts, err)

    plt.show()
    
    print(time.time() - strart_time)
    pass
