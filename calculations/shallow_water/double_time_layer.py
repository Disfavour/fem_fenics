from fenics import *
import numpy as np
from scipy.constants import g
import analytic


set_log_level(LogLevel.WARNING)


def calculate(hl, hr, T, mesh_size, tau, theta, sigma, ts2store, info=False):
    domain_size = 5
    degree = 1
    mesh_size_exact = 100000

    m, E, m_e, E_e  = [], [], [], []
    h_, u_, h_e, u_e = [], [], [], []
    err_h, err_u = [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()
    mesh_exact = IntervalMesh(mesh_size_exact, -domain_size, domain_size)
    x_e = mesh_exact.coordinates().flatten()

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))
    H1 = FiniteElement("P", mesh.ufl_cell(), 1)
    U1 = FiniteElement("P", mesh.ufl_cell(), 1)
    WE = FunctionSpace(mesh_exact, MixedElement([H1, U1]))

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)
    we = Function(WE)
    he, ue = split(we)

    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx
    m_eq_exact = he * dx
    E_eq_exact = 0.5*(he*ue*ue + g*he*he) * dx

    h = np.array((h, hn))
    u = np.array((u, un))

    dt = lambda f: (f[0]-f[1])/tau
    s = lambda f: sigma*f[0] + (1-sigma)*f[1]
    
    F = (dt(h) + s(h*u).dx(0)) * ht*dx \
        + (dt(h*u) + s(h*u**2).dx(0) + g/2*s(h**2).dx(0)) * ut*dx

    get_exact_solution = analytic.get_solution_func(hl, hr)

    def set_exact(f, x, FS):
        hs_exact, us_exact = get_exact_solution(x, t)
        f.vector()[dof_to_vertex_map(FS)[FS.sub(0).dofmap().dofs()]] = hs_exact
        f.vector()[dof_to_vertex_map(FS)[FS.sub(1).dofmap().dofs()]] = us_exact

    def collect_data():
        set_exact(we, x_e, WE)

        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        m_e.append(assemble(m_eq_exact))
        E_e.append(assemble(E_eq_exact))
        err_h.append(errornorm(w.sub(0), we.sub(0), norm_type='L2', degree_rise=0))
        err_u.append(errornorm(w.sub(1), we.sub(1), norm_type='L2', degree_rise=0))

        if np.isclose(t, ts2store).any():
            h_.append(w.sub(0).compute_vertex_values())
            u_.append(w.sub(1).compute_vertex_values())
            h_e.append(we.sub(0).compute_vertex_values())
            u_e.append(we.sub(1).compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

    t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=1), W))
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

    return map(np.array, (ts, m, E, m_e, E_e, err_h, err_u, x, h_, u_, x_e, h_e, u_e))


if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    strart_time = time.time()
    #ts, m, E, m_e, E_e, err_h, err_u, x, h_, u_, x_e, h_e, u_e = calculate(2, 1, 0.9, 400, 0.01, None, 1.0, [0.9], True)
    ts, m, E, m_e, E_e, err_h, err_u, x, h_, u_, x_e, h_e, u_e = calculate(2, 1, 0.9, 200, 0.005, None, 0.9, [0.3, 0.6, 0.9], True)
    plt.plot(x, h_.T)
    #plt.plot(x_e, h_e[0], ':')
    plt.show()
    print(time.time() - strart_time)
    pass
