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
        + kappa * (u-un).dx(0) * ut*dx(mesh)

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


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    strart_time = time.time()
    # mesh_size=200, tau=0.0025
    # mesh_size=100, tau=0.005

    # errs = []

    # for k in np.arange(0, 0.08, 0.01):
    #     print(k, end=' ')
    #     x, u, t, err, x_e, u_exact = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=k, ntest=2, info=False)
    #     errs.append(max(err))
    #     print(errs[-1])

    # for i, e in enumerate(errs):
    #     print(i, e)
    # #print(errs)
    # exit()

    x, u1, t, err1, x_e, u_exact = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.0, ntest=1, info=True)
    x, u2, t, err2, x_e, u_exact = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.04, ntest=1, info=True)
    x, u3, t, err3, x_e, u_exact = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.08, ntest=1, info=True)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    for u, c in zip([u1, u2, u3], ('r', 'g', 'b')):
        for uc in u:
            plt.plot(x, uc, c)

    # for un1, un2, un3, u_e, c in zip(u1, u2, u3, u_exact, ('r', 'g', 'b', 'c')):
    #     plt.plot(x, un1, '-.' + c)
    #     plt.plot(x, un2, '-' + c)
    #     plt.plot(x, un3, '--' + c)
        #plt.plot(x_e, u_e, ':' + c)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(t, err1)
    plt.plot(t, err2)
    plt.plot(t, err3)
    plt.legend(["1", '2', '3'])

    plt.show()
    
    # x1, u1, t2, err2 = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.08, ntest=1, info=True)

    # x2, u2, t3, err3 = calculate(mesh_size=200, tau=0.0025, degree=2, T=0.6, ts2store=[0.0, 0.2, 0.4, 0.6], kappa=0.16, ntest=1, info=True)

    # plt.figure()
    # plt.plot(x, u.T, '-r')
    # plt.plot(x1, u1.T, '--b')
    # plt.plot(x2, u2.T, ':g')

    # plt.figure()
    # plt.plot(t1, err1)
    # plt.plot(t2, err2)
    # plt.plot(t3, err3)
    # plt.legend(["1", '2', '3'])

    # plt.figure()
    # plt.plot(t, m)
    # plt.ylim(0, 2)
    #plt.savefig('1234')
    #plt.show()

    # e^(-100800*(x-0.1)^(4))
    
    print(time.time() - strart_time)
    pass
