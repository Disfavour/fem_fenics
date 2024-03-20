from fenics import *
import numpy as np
from scipy.constants import g
import analytic


set_log_level(LogLevel.WARNING)


def calculate(mesh_size, tau, degree, T, ts2store, info=False):
    domain_size = 5

    u_= []

    ts = np.arange(0, T+tau/2, tau)
    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)

    bc = DirichletBC(P, Constant(0), CompiledSubDomain('on_boundary && near(x[0], left)', left=-domain_size))

    un = Function(P)
    u = TrialFunction(P)
    ut = TestFunction(P)

    F = ((u - un)/tau + u.dx(0)) * ut*dx
    a, L = lhs(F), rhs(F)

    u = Function(P)

    def collect_data():
        if np.isclose(t, ts2store).any():
            u_.append(u.compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f}')

    t = 0
    u.assign(project(Expression('x[0] <= -4 ? 1 : 0', degree=1), P))
    collect_data()
    un.assign(u)

    for t in ts[1:]:
        solve(a == L, u, bc)
        collect_data()
        un.assign(u)

    return map(np.array, (x, u_))


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    strart_time = time.time()
    x, u = calculate(mesh_size=400, tau=0.01, degree=1, T=1.0, ts2store=[0.3, 0.5, 0.9], info=True)
    plt.figure()
    plt.plot(x, u.T)
    plt.figure()
    for i in u:
        plt.plot(x, i)
    plt.show()
    #calculate(2, 1, 0.9, 400, 0.01, 0.1, 1.0, [], True)
    print(time.time() - strart_time)
    pass
