from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def shallow_water_1d(s=1, tau=0.01, mesh_size=100, T=1, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None

    hl, hr = 10, 0.7
    domain_size = 5

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    U = FiniteElement('P', mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    w.rename('(h, u)', 'shallow water 1d')

    def hs():
        return s*h + (1-s)*hn

    def us():
        return s*u + (1-s)*un
    
    # F = ((h - hn)/tau + div(h*u)) * ht * dx \
    #      + ((h*u - hn*un)/tau + div(h*outer(u, u)) + g*h*grad(h)) * ut * dx

    # F = ((h-hn)/tau + (hs()*us()).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs()*us()*us()).dx(0) + g/2*(hs()*hs()).dx(0)) * ut*dx

    F = (h - hn)/tau*ht*dx - hs()*us()*ht.dx(0)*dx \
        + (h*u - hn*un)/tau*ut*dx - hs()*us()*us()*ut.dx(0)*dx + g*hs()*hs().dx(0)*ut*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        if vtkfile is not None:
            vtkfile << (w, t)

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')
        print(w.sub(0).compute_vertex_values().min())

    # t = 0
    w.assign(project(Expression(('x[0] < 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bc)

        collect_data()

        wn.assign(w)


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(dirname(__file__))
    paraview = join(base_dir, 'paraview')
    vtkfname = join(paraview, 'numerical.pvd')
    shallow_water_1d(s=1.0, tau=0.01, mesh_size=200, vtkfname=None)
