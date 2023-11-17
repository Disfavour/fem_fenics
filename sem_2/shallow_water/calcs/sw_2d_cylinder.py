from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

#set_log_level(LogLevel.WARNING)


def shallow_water_1d(s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_size = 5

    vtkfile = File(join(paraview, 'sw_2d_cylinder.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    U = VectorElement('P', mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bcs = [
        DirichletBC(W.sub(1).sub(0), Constant(0), CompiledSubDomain("on_boundary && (near(x[0], -xsize) || near(x[0], xsize))", xsize=domain_size)),
        DirichletBC(W.sub(1).sub(1), Constant(0), CompiledSubDomain("on_boundary && (near(x[1], -ysize) || near(x[1], ysize))", ysize=domain_size)),
    ]

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    w.rename('(h, u)', 'shallow water 2d cylinder')

    def hs():
        return s*h + (1-s)*hn

    def us():
        return s*u + (1-s)*un
    
    F = (h - hn)/tau*ht*dx - hs()*dot(us(), grad(ht))*dx \
        + dot((h*u - hn*un)/tau, ut)*dx \
        - inner(hs()*outer(us(), us()), nabla_grad(ut)) * dx \
        + dot(g*hs()*grad(h), ut)*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0]*x[0] + x[1]*x[1] < 2.5 ? hl : hr', '0', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bcs)

        collect_data()

        wn.assign(w)


if __name__ == '__main__':
    shallow_water_1d(s=0.5, tau=0.01, mesh_size=100)
