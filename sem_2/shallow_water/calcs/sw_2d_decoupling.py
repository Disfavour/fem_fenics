from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

set_log_level(LogLevel.WARNING)


def shallow_water_1d(K=1, s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_x = 5
    domain_y = 2
    ms_y = round(domain_y / domain_x * mesh_size)

    vtkfile = File(join(paraview, 'sw_2d_decoupling.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = RectangleMesh(Point(-domain_x, -domain_y), Point(domain_x, domain_y), mesh_size, ms_y)

    H = FunctionSpace(mesh, 'P', degree)
    U = VectorFunctionSpace(mesh, 'P', degree)

    #W = VectorFunctionSpace(mesh, 'P', degree, 3)
    H_ = FiniteElement('P', mesh.ufl_cell(), degree)
    U_ = VectorElement('P', mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H_, U_]))

    bcs = [
        DirichletBC(U.sub(0), Constant(0), CompiledSubDomain("on_boundary && (near(x[0], -xsize) || near(x[0], xsize))", xsize=domain_x)),
        DirichletBC(U.sub(1), Constant(0), CompiledSubDomain("on_boundary && (near(x[1], -ysize) || near(x[1], ysize))", ysize=domain_y)),
    ]

    ht, ut = TestFunction(H), TestFunction(U)
    h, u_ = TrialFunction(H), TrialFunction(U)
    u = Function(U)
    hn, un = Function(H), Function(U)

    w = Function(W)
    w.rename('(h, u)', 'shallow water 2d decoupling')

    def hs():
        return s*h + (1-s)*hn
    
    F1 = (h - hn)/tau*ht*dx - hs()*dot(u, grad(ht))*dx
    a1, L1 = lhs(F1), rhs(F1)

    h = Function(H)

    def us():
        return s*u_ + (1-s)*un

    F2 = dot((h*u_ - hn*un)/tau, ut)*dx - inner(h*outer(us(), u), nabla_grad(ut))*dx + dot(g*h*grad(h), ut)*dx
    a2, L2 = lhs(F2), rhs(F2)
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        assign(w.sub(0), h)
        assign(w.sub(1), u)

        vtkfile << (w, t)

    # t = 0
    h.assign(project(Expression('x[0] < 0 ? hl : hr', hl=hl, hr=hr, degree=degree), H))
    u.assign(project(Constant((0, 0)), U))

    collect_data()

    hn.assign(h)
    un.assign(u)

    for t in ts[1:]:
        for k in range(K):
            solve(a1 == L1, h)
            solve(a2 == L2, u, bcs)

        collect_data()

        hn.assign(h)
        un.assign(u)


if __name__ == '__main__':
    shallow_water_1d(K=1, s=1, tau=0.01, mesh_size=100)
