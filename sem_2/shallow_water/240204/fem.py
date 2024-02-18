from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def get_solution(hl, hr, T, mesh_size, tau, s, ts2store):
    domain_size = 5
    degree = 1
    ms, Es = [], []
    res_h, res_u = [], []

    t = 0
    ts = np.arange(t, T+tau/2, tau)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)
    x = mesh.coordinates().flatten()

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    # F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * ut*dx

    F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g*hs*(hs).dx(0)) * ut*dx
    
    # F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0)) * ut*dx - g/2*(hs*hs) * ut.dx(0) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))
        ms.append(m)
        Es.append(E)
        print(f'Time {t:>7.5f} m {m:>7.5f} E {E:>7.5f}')

        if np.isclose(t, ts2store).any():
            res_h.append(w.sub(0).compute_vertex_values())
            res_u.append(w.sub(1).compute_vertex_values())

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

    return ts, ms, Es, x, res_h, res_u
    #return map(np.array, (ts, ms, Es, x, res_h))


if __name__ == '__main__':
    get_solution(hl=2, hr=1, T=0.5, mesh_size=200, tau=0.005, s=1.0, ts2store=[])
    pass
