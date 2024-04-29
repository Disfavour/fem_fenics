from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def calculate(hl, hr, T, M, tau, fname, info=False):
    vtkfile = File(fname)
    degree = 1

    x0, y0 = -5, -5
    x1, y1 = 5, 5

    m, E = [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), M, M)

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = VectorElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))
    
    bc = [
        DirichletBC(W.sub(1).sub(0), Constant(0), CompiledSubDomain('on_boundary && (near(x[0], left) || near(x[0], right))', left=x0, right=x1)),
        DirichletBC(W.sub(1).sub(1), Constant(0), CompiledSubDomain('on_boundary && (near(x[1], top) || near(x[1], bot))', top=y1, bot=y0)),
    ]

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx

    F = (h-hn)/tau * ht*dx - dot(h*u, grad(ht)) * dx \
        + dot((h*u - hn*un)/tau, ut) * dx \
        - inner(outer(h*u, u), nabla_grad(ut)) * dx \
        + g * dot(h * grad(h), ut) * dx
    
    def collect_data():
        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        vtkfile << (w, t)

        if info:
            print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

    t = 0
    w.assign(project(Expression(('1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))', '0', '0'), alf=2, bet=20, degree=degree), W))
    collect_data()
    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bc, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-13,
            'relative_tolerance': 1e-13,
            'maximum_iterations': 50,
            'relaxation_parameter': 1.0,
        }})
        collect_data()
        wn.assign(w)

    return map(np.array, (ts, m, E))


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    from os.path import join

    strart_time = time.time()
    t, m, E = calculate(2, 1, 2.0, 200, 0.005, join('data', 'paraview', 'sw_2d_2_.pvd'), True)
    print(time.time() - strart_time)

    plt.figure()
    plt.plot(t, m)

    plt.figure()
    plt.plot(t, E)

    plt.show()
