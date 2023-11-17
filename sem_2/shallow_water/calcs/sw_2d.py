from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

#set_log_level(LogLevel.WARNING)
set_log_level(LogLevel.PROGRESS)


def shallow_water_1d(s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_x = 5
    domain_y = 1
    ms_y = round(domain_y / domain_x * mesh_size)

    vtkfile = File(join(paraview, f'plot_sw_2d_s{s}_tau{tau}_ms{mesh_size}_.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = RectangleMesh(Point(-domain_x, -domain_y), Point(domain_x, domain_y), mesh_size, ms_y)

    H = FiniteElement('P', mesh.ufl_cell(), degree)
    U = VectorElement('P', mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bcs = [
        DirichletBC(W.sub(1).sub(0), Constant(0), CompiledSubDomain("on_boundary && (near(x[0], -xsize) || near(x[0], xsize))", xsize=domain_x)),
        DirichletBC(W.sub(1).sub(1), Constant(0), CompiledSubDomain("on_boundary && (near(x[1], -ysize) || near(x[1], ysize))", ysize=domain_y)),
    ]

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    w.rename('(h, u)', 'shallow water 2d')

    def hs():
        return s*h + (1-s)*hn

    def us():
        return s*u + (1-s)*un
    
    F = (h - hn)/tau*ht*dx - hs()*dot(us(), grad(ht))*dx \
        + dot((h*u - hn*un)/tau, ut)*dx \
        - inner(hs()*outer(us(), us()), nabla_grad(ut))*dx \
        + dot(g*hs()*grad(hs()), ut)*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0] < 0 ? hl : hr', '0', '0'), hl=hl, hr=hr, degree=0), W))

    collect_data()

    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bcs)

        collect_data()

        wn.assign(w)


if __name__ == '__main__':
    shallow_water_1d(s=1, tau=0.005, mesh_size=200)
    #shallow_water_1d(s=1, tau=0.005, mesh_size=200)

    # from multiprocessing import Pool

    # # p1 = list((1, tau, 200) for tau in (0.01, 0.005, 0.0025))
    # # p2 = list((s, 0.005, 200) for s in (0.75, 1.5))
    # # p3 = list((1, 0.005, ms) for ms in (100, 400))

    # # t / 2
    # p1 = list((1, tau, 200) for tau in (0.005, 0.0025, 0.00125))
    # p2 = list((s, 0.0025, 200) for s in (0.75, 1.5))
    # p3 = list((1, 0.0025, ms) for ms in (100, 400))

    # params = p1 + p2 + p3
    # # print(params)
    # # exit()

    # #params = ((s, tau, ms) for ms in [400, 200, 100] for tau in [0.0025, 0.005, 0.01] for s in [0.75, 1, 1.5]) 

    # with Pool() as p:
    #     res = p.starmap_async(shallow_water_1d, params)
    #     res.wait()
