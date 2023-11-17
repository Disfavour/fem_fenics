from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

set_log_level(LogLevel.WARNING)
#set_log_level(LogLevel.PROGRESS)


def sw_2d_w(s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_x = 5
    domain_y = 1
    ms_y = round(domain_y / domain_x * mesh_size)

    vtkfile = File(join(paraview, f'plot_sw_2d_w_s{s}_tau{tau}_ms{mesh_size}_.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = RectangleMesh(Point(-domain_x, -domain_y), Point(domain_x, domain_y), mesh_size, ms_y)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    W = VectorElement('P', mesh.ufl_cell(), degree)
    X = FunctionSpace(mesh, MixedElement([H, W]))

    bcs = [
        DirichletBC(X.sub(1).sub(0), Constant(0), CompiledSubDomain("on_boundary && (near(x[0], -xsize) || near(x[0], xsize))", xsize=domain_x)),
        DirichletBC(X.sub(1).sub(1), Constant(0), CompiledSubDomain("on_boundary && (near(x[1], -ysize) || near(x[1], ysize))", ysize=domain_y)),
    ]

    # w = h*u
    x, xn = Function(X), Function(X)
    ht, wt = TestFunctions(X)
    h, w = split(x)
    hn, wn = split(xn)

    x.rename(r'(h, h*u)', 'shallow water 2d')

    def hs():
        return s*h + (1-s)*hn

    def ws():
        return s*w + (1-s)*wn
    
    F = (h - hn)/tau*ht*dx - dot(ws(), grad(ht))*dx \
        + dot((w - wn)/tau, wt)*dx \
        - inner(outer(ws(), ws() / hs()), nabla_grad(wt)) * dx \
        + dot(g*hs()*grad(hs()), wt)*dx
    
    # F = (h - hn)/tau*ht*dx - hs()*dot(us(), grad(ht))*dx \ #+ div(ws())*ht*dx \
    #     + dot((h*u - hn*un)/tau, ut)*dx \
    #     - inner(hs()*outer(us(), us()), nabla_grad(ut)) * dx \
    #     + dot(g*hs()*grad(h), ut)*dx
    
    m_eq = h * dx
    E_eq = 0.5*(dot(w, w) / h + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        vtkfile << (x, t)

    # t = 0
    x.assign(project(Expression(('x[0] < 0 ? hl : hr', '0', '0'), hl=hl, hr=hr, degree=degree), X))

    collect_data()

    xn.assign(x)

    for t in ts[1:]:
        solve(F == 0, x, bcs)

        collect_data()

        xn.assign(x)


if __name__ == '__main__':
    sw_2d_w(s=0.5, tau=0.01, mesh_size=100)
    # sw_2d_w(s=1, tau=0.005, mesh_size=400)
