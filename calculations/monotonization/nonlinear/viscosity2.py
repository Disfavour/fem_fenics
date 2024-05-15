from fenics import *
import numpy as np


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, kappa, ts2store, ntest=2, info=False):
    sigma = 0.5
    degree = 1
    T = 0.6
    u_, err = [], []
    ue_ = []

    if ntest == 1:
        a, b, x0 = 0.5, 1.5, 0.2
    elif ntest == 2:
        a, b, x0 = 1.5, 0.5, 0.2

    ts = np.arange(0, T+tau/2, tau)
    mesh = UnitIntervalMesh(mesh_size)
    x = mesh.coordinates().flatten()
    mesh_e = UnitIntervalMesh(10000)
    xe = mesh_e.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)
    PE = FunctionSpace(mesh_e, 'P', degree)

    bc = DirichletBC(P, Constant(a), CompiledSubDomain('on_boundary && near(x[0], 0)'))
    
    u, un = Function(P), Function(P)
    ut = TestFunction(P)
    u_e = Function(PE)

    dt = (u - un)/tau
    us = sigma*u + (1-sigma)*un

    # F = (dt + us*us.dx(0)) * ut*dx \
    #     + kappa * us**2 * (u-un).dx(0)*ut.dx(0) * dx(mesh)
    
    # справедливо только для линейных операторов, поэтому us**2 убираем
    F = (dt + us*us.dx(0)) * ut*dx \
        + kappa * (u-un).dx(0)*ut.dx(0) * dx(mesh)
    
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
        u_e.vector()[dof_to_vertex_map(PE)[PE.dofmap().dofs()]] = exact(t)
        err.append(errornorm(u, u_e, norm_type='L2', degree_rise=0))

        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())
            ue_.append(u_e.compute_vertex_values())
        
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
    
    x, u, x_e, ue, t, err = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.2, 0.4, 0.6], kappa=0.02, ntest=2, info=True)
    x, u2, x_e, ue, t, err2 = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.2, 0.4, 0.6], kappa=0.04, ntest=2, info=True)
    x, u3, x_e, ue, t, err3 = calculate(mesh_size=200, tau=0.0025, degree=1, T=0.6, ts2store=[0.2, 0.4, 0.6], kappa=0.08, ntest=2, info=True)

    lines = ['-', '--', '-.', ':']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    for u_e, u1_, u2_, u3_, c in zip(ue, u, u2, u3, colors):
        plt.plot(x_e, u_e, c)
        plt.plot(x, u1_, '--'+c)
        plt.plot(x, u2_, '-.'+c)
        plt.plot(x, u3_, ':'+c)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
    plt.plot(t, err)
    plt.plot(t, err2)
    plt.plot(t, err3)
    plt.legend(['1', '2', '3'])

    plt.show()
    
    print(time.time() - strart_time)
    pass
