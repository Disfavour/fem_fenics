from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def sw_1d_smooth_v2(hl, s=1, tau=0.01, mesh_size=100, T=1, dif_t=[], degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None
    domain_size = 5
    Es = []
    res_h, res_u = [], []
    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    U = FiniteElement('P', mesh.ufl_cell(), degree)
    Q = FunctionSpace(mesh, MixedElement([H, U]))

    bc = DirichletBC(Q.sub(1), Constant(0), 'on_boundary')

    q, qn = Function(Q), Function(Q)
    ht, ut = TestFunctions(Q)
    h, u = split(q)
    hn, un = split(qn)

    q.rename('(h, u)', 'shallow water 1d')

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un
    
    F = (h - hn)/tau*ht*dx + (hs*us).dx(0)*ht*dx \
        + (h*u - hn*un)/tau*ut*dx \
        + (hs*us*us).dx(0)*ut*dx \
        + g*hs*hs.dx(0)*ut*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))
        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        Es.append(E)

        if np.isclose(t, dif_t).any():
            res_h.append(q.sub(0).compute_vertex_values())
            res_u.append(q.sub(1).compute_vertex_values())
        
        if vtkfile is not None:
            vtkfile << (q, t)

    # t = 0
    q.assign(project(Expression(('x[0] <= 0 ? hl : 1 + alf*exp(-bet*(x[0]*x[0]))', '0'), hl=hl, alf=hl-1, bet=20, degree=degree), Q))
    collect_data()
    qn.assign(q)

    for t in ts[1:]:
        solve(F == 0, q, bc)
        collect_data()
        qn.assign(q)
    
    return map(np.array, (mesh.coordinates(), res_h, res_u, ts, Es))


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(__file__)
    paraview = join(base_dir, 'paraview')
    vtkfname = join(paraview, 'sw1d_s1_v1_0.055.pvd')
    sw_1d_smooth_v2(hl=3, vtkfname=None)
