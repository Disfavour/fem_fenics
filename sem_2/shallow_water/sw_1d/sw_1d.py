from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def sw_1d(s=1, tau=0.01, mesh_size=100, T=1, time_moments=[], depth=1, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None
    domain_size = 5
    Es = []
    res_h, res_u = [], []
    courant = []
    delta = [0]
    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    delta_x = mesh.coordinates().flatten()[1] - mesh.coordinates().flatten()[0]

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

    # F = (h - hn)/tau*ht*dx - hs*us*ht.dx(0)*dx \
    #     + (h*u - hn*un)/tau*ut*dx - hs*us*us*ut.dx(0)*dx + g*hs*hs.dx(0)*ut*dx
    
    F = (h - hn)/tau*ht*dx + (hs*us).dx(0)*ht*dx \
        + (h*u - hn*un)/tau*ut*dx \
        + (hs*us*us).dx(0)*ut*dx \
        + g*hs*hs.dx(0)*ut*dx
    
    #hts = s*ht + (1-s)*hn

    #P = FunctionSpace(mesh, 'P', 1)

    #delta_c = Function(Q)
    #dt, d_t = TestFunctions(Q)
    #d, delta_ = split(delta_c)
    
    #F_d = ((h - hn)/tau*us*us*dx + (hs*us).dx(0)*us*us*dx) + 2*d
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))
        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        #solve(delta==0, )
        delta.append(delta[-1] + tau*assemble(-0.5*((h - hn)/tau*us*us*dx + (hs*us).dx(0)*us*us*dx)))

        Es.append(E)

        courant.append(abs(q.sub(1).compute_vertex_values()).max() * tau / delta_x)

        if np.isclose(t, time_moments).any():
            res_h.append(q.sub(0).compute_vertex_values())
            res_u.append(q.sub(1).compute_vertex_values())
        
        if vtkfile is not None:
            vtkfile << (q, t)
        
        if (q.sub(0).compute_vertex_values().min() < 0):
            print(q.sub(0).compute_vertex_values().min())

    # t = 0
    q.assign(project(Expression(('depth + alf*exp(-bet*(x[0]*x[0]))', '0'), depth=depth, alf=2, bet=20, degree=degree), Q))
    collect_data()
    qn.assign(q)

    for t in ts[1:]:
        solve(F == 0, q, bc)
        collect_data()
        qn.assign(q)
    
    return map(np.array, (mesh.coordinates(), res_h, res_u, time_moments, ts, Es, courant, delta))


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(__file__)
    paraview = join(base_dir, 'paraview')
    vtkfname = join(paraview, 'sw1d_s1_v1_0.055.pvd')
    sw_1d(s=1.0, tau=0.01, mesh_size=200, T=5, depth=1, vtkfname=None)
