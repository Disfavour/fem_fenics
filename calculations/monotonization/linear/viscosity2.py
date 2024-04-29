from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, degree, T, ts2store, kappa, ntest, info=False):
    sigma = 0.5
    u_, err = [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = UnitIntervalMesh(mesh_size)
    x = mesh.coordinates().flatten()
    mesh_e = UnitIntervalMesh(10000)

    P = FunctionSpace(mesh, 'P', degree)
    PE = FunctionSpace(mesh_e, 'P', degree)

    bc = DirichletBC(P, Constant(0), CompiledSubDomain('on_boundary && near(x[0], 0)'))
    
    u, un = Function(P), Function(P)
    ut = TestFunction(P)
    u_e = Function(PE)

    dt = (u - un)/tau
    us = sigma*u + (1-sigma)*un

    F = (dt + us.dx(0)) * ut*dx \
        + kappa * (u-un).dx(0)*ut.dx(0) * dx(mesh)

    def collect_data():
        u_expr.t = t
        u_e.assign(project(u_expr, PE))
        err.append(errornorm(u, u_e, norm_type='L2', degree_rise=0))

        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())

        if info:
            print(f'Time {t:>7.5f}')

    t = 0

    if ntest == 1:
        u_expr = Expression('x[0] >= 0.15 + t && x[0] <= 0.25 + t ? 1 : 0', t=t, degree=1)
    elif ntest == 2:
        u_expr = Expression('x[0] >= 0.15 + t && x[0] <= 0.25 + t ? sin((x[0] - 0.15 - t) * 10 * pi) : 0', t=t, degree=degree)
    elif ntest == 3:
        u_expr = Expression('exp(-1000*pow(x[0] - 0.2 - t, 2))', t=t, degree=degree)
    
    u.assign(project(u_expr, P))
        
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

    return map(np.array, (x, u_, ts, err))

def calculate(mesh_size, tau, degree, T, ts2store, kappa, ntest, info=False):
    sigma = 0.5
    u_, m = [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = UnitIntervalMesh(mesh_size)
    x = mesh.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)

    bc = DirichletBC(P, Constant(0), CompiledSubDomain('on_boundary && near(x[0], 0)'))
    
    u, un = Function(P), Function(P)
    ut = TestFunction(P)

    dt = (u - un)/tau
    us = sigma*u + (1-sigma)*un

    F = (dt + us.dx(0)) * ut*dx \
        + kappa * (u-un).dx(0)*ut.dx(0) * dx(mesh)

    def collect_data():
        m.append(assemble(u * dx))
        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f} m {m[-1]:>7.5f}')

    t = 0

    if ntest == 1:
        u.assign(project(Expression('x[0] >= 0.15 && x[0] <= 0.25 ? 1 : 0', degree=1), P))
    elif ntest == 2:
        u.assign(project(Expression('x[0] >= 0.15 && x[0] <= 0.25 ? sin((x[0] - 0.15) * 10 * pi) : 0', degree=degree), P))
    elif ntest == 3:
        u.assign(project(Expression('exp(-2000*pow(x[0] - 0.2, 2))', degree=degree), P))
        
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
    # mesh_size=200, tau=0.0025
    # mesh_size=100, tau=0.005

    x, u, t, m = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.0, ntest=1, info=True)
    
    x1, u1, t, m = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.04, ntest=1, info=True)

    x2, u2, t, m = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.08, ntest=1, info=True)

    plt.figure()
    plt.plot(x, u.T, '-r')
    plt.plot(x1, u1.T, '--b')
    plt.plot(x2, u2.T, ':g')

    # plt.figure()
    # plt.plot(t, m)
    # plt.ylim(0, 2)
    #plt.savefig('1234')
    plt.show()

    # e^(-100800*(x-0.1)^(4))
    
    print(time.time() - strart_time)
    pass
