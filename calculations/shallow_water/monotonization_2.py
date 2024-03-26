from fenics import *
import numpy as np
from scipy.constants import g
import analytic


set_log_level(LogLevel.WARNING)


def calculate(hl, hr, T, mesh_size, tau, alfa, gamma, ts2store, info=False):
    domain_size = 5
    degree = 1
    mesh_size_exact = 10000

    m, E, m_e, E_e  = [], [], [], []
    h_, u_, h_e, u_e = [], [], [], []
    err_h, err_u = [], []

    ts = np.arange(0, T+tau/2, tau)
    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()
    mesh_exact = IntervalMesh(mesh_size_exact, -domain_size, domain_size)
    x_e = mesh_exact.coordinates().flatten()

    P = FunctionSpace(mesh, 'P', degree)
    PE = FunctionSpace(mesh_exact, 'P', degree)

    h, hn, u = Function(P), Function(P), Function(P)
    ht = TestFunction(P)
    he, ue = Function(PE), Function(PE)

    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx
    m_eq_exact = he * dx
    E_eq_exact = 0.5*(he*ue*ue + g*he*he) * dx

    dth = (h - hn) / tau
    b = abs(h.dx(0)) ** gamma
    
    # F = (dth + (h*u).dx(0)) * ht*dx 
    # F = (dth - a*tau**2*(b*h.dx(0)).dx(0) + (h*u).dx(0)) * ht*dx
    # F = (dth + (h*u).dx(0)) * ht*dx - a*tau**2*(b*ht*h.dx(0)).dx(0)*dx + a*tau**2*b*h.dx(0)*ht.dx(0)*dx
    F = (dth + (h*u).dx(0)) * ht*dx + alfa*tau**2*b*h.dx(0)*ht.dx(0)*dx(mesh)

    get_exact_solution = analytic.get_solution_func(hl, hr)

    def collect_data():
        h_exact, u_exact = get_exact_solution(x_e, t)
        he.vector()[dof_to_vertex_map(PE)] = h_exact
        ue.vector()[dof_to_vertex_map(PE)] = u_exact

        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        m_e.append(assemble(m_eq_exact))
        E_e.append(assemble(E_eq_exact))
        err_h.append(errornorm(h, he, norm_type='L2', degree_rise=0))
        err_u.append(errornorm(u, ue, norm_type='L2', degree_rise=0))

        if np.isclose(t, ts2store).any():
            h_.append(h.compute_vertex_values())
            u_.append(u.compute_vertex_values())
            h_e.append(he.compute_vertex_values())
            u_e.append(ue.compute_vertex_values())
        
        if info:
            print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

    t = 0
    h.assign(project(Expression('x[0] <= 0 ? hl : hr', hl=hl, hr=hr, degree=1), P))
    u.assign(project(Constant(0), P))
    collect_data()
    hn.assign(h)

    for t in ts[1:]:
        solve(F == 0, h, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-13,
            'relative_tolerance': 1e-13,
            'maximum_iterations': 50,
            'relaxation_parameter': 1.0,
        }})
        u.vector()[dof_to_vertex_map(P)] = get_exact_solution(x, t)[1]
        collect_data()
        hn.assign(h)

    return map(np.array, (ts, m, E, m_e, E_e, err_h, err_u, x, h_, u_, x_e, h_e, u_e))


if __name__ == '__main__':
    import time
    strart_time = time.time()
    calculate(2, 1, 0.9, 400, 0.01, 0.1, 1.0, [], True)
    print(time.time() - strart_time)
    pass
