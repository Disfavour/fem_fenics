from fenics import *
import numpy as np
from scipy.constants import g
import analytic


set_log_level(LogLevel.WARNING)


def get_solution(hl, hr, T, mesh_size, tau, s, ts2store):
    domain_size = 5
    degree = 1
    mesh_size_exact = 50000

    m, E = [], []
    h_, u_ = [], []
    m_e, E_e = [], []
    h_e, u_e = [], []
    err_h, err_u = [], []

    t = 0
    ts = np.arange(t, T+tau/2, tau)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()

    mesh_exact = IntervalMesh(mesh_size_exact, -domain_size, domain_size)
    x_e = mesh_exact.coordinates().flatten()

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    We = FunctionSpace(mesh_exact, MixedElement([H, U]))
    We_component_idxs = (
        We.sub(0).dofmap().dofs(),
        We.sub(1).dofmap().dofs())
    We_d2v_for_components = (
        dof_to_vertex_map(We)[We_component_idxs[0]],
        dof_to_vertex_map(We)[We_component_idxs[1]])

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    we = Function(We)
    he, ue = split(we)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g*hs*(hs).dx(0)) * ut*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    m_eq_exact = he * dx
    E_eq_exact = 0.5*(he*ue*ue + g*he*he) * dx

    get_exact_solution = analytic.get_solution_func(hl, hr)

    def collect_data():
        hs_exact, us_exact = get_exact_solution(x_e, t)
        we.vector()[We_d2v_for_components[0]] = hs_exact
        we.vector()[We_d2v_for_components[1]] = us_exact

        m.append(assemble(m_eq))
        E.append(assemble(E_eq))

        m_e.append(assemble(m_eq_exact))
        E_e.append(assemble(E_eq_exact))

        print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

        if np.isclose(t, ts2store).any():
            h_.append(w.sub(0).compute_vertex_values())
            u_.append(w.sub(1).compute_vertex_values())

            h_e.append(we.sub(0).compute_vertex_values())
            u_e.append(we.sub(1).compute_vertex_values())
        
        err_h.append(errornorm(we.sub(0), w.sub(0), norm_type='L2', degree_rise=0, mesh=mesh_exact))
        err_u.append(errornorm(we.sub(1), w.sub(1), norm_type='L2', degree_rise=0, mesh=mesh_exact))

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=1), W))
    collect_data()
    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bc, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-8,
            'relative_tolerance': 1e-8,
            'maximum_iterations': 50,
            'relaxation_parameter': 1.0,
        }})
        collect_data()
        wn.assign(w)

    return map(np.array, (ts, m, E, m_e, E_e, err_h, err_u, x, h_, u_, x_e, h_e, u_e))


if __name__ == '__main__':
    import time
    strart_time = time.time()
    get_solution(10, 1, 0.5, 200, 0.005, 1, [])
    print(time.time() - strart_time)
    pass
