from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

set_log_level(LogLevel.WARNING)


def shallow_water_1d(K=5, s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_size = 5

    #vtkfile = File(join(paraview, 'decoupling.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FunctionSpace(mesh, 'P', degree)
    U = FunctionSpace(mesh, 'P', degree)

    bc = DirichletBC(U, Constant(0), 'on_boundary')

    ht, ut = TestFunction(H), TestFunction(U)
    h, u = TrialFunction(H), TrialFunction(U)
    hn, un = Function(H), Function(U)
    hk, uk = Function(H), Function(U)

    #h.rename('(h, u)', 'shallow water 1d')

    # def hs():
    #     return s*h_ + (1-s)*hn
    
    # def hsk():
    #     return s*h + (1-s)*hn

    # def us():
    #     return s*u_ + (1-s)*un
    
    # def usk():
    #     return s*uk + (1-s)*un
    #parameters["linear_algebra_backend"] = "PETSc"

    #F1 = ((h - hn)/tau + (h*uk).dx(0)) * ht*dx
    F1 = (h - hn)/tau*ht*dx - h*uk*ht.dx(0)*dx
    a1, L1 = lhs(F1), rhs(F1)

    h = Function(H)

    #F2 = ((h*u - hn*un)/tau + (h*uk*u).dx(0) + g*h*h.dx(0)) * ut*dx
    F2 = (h*u - hn*un)/tau*ut*dx - h*uk*u*ut.dx(0)*dx + g*h*h.dx(0)*ut*dx
    a2, L2 = lhs(F2), rhs(F2)

    u = Function(U)
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.15f} E {E:>7.15f}')

        #vtkfile << (h, t)

    # t = 0
    h.assign(project(Expression('x[0] < 0 ? hl : hr', hl=hl, hr=hr, degree=degree), H))
    u.assign(project(Constant(0), U))

    collect_data()

    hn.assign(h)
    un.assign(u)

    for t in ts[1:]:
        hk.assign(hn)
        uk.assign(un)
        
        for k in range(K):
            # # Assemble system
            # A = assemble(a1)
            # b = assemble(L1)
            # solver = PETScKrylovSolver("cg")
            # solver.set_operator(A)
            # null_vec = Vector(h.vector())
            # H.dofmap().set(null_vec, 1.0)
            # null_vec *= 1.0/null_vec.norm("l2")
            # null_space = VectorSpaceBasis([null_vec])
            # as_backend_type(A).set_nullspace(null_space)
            # null_space.orthogonalize(b)
            # solver.solve(h.vector(), b)

            solve(a1 == L1, h)
            solve(a2 == L2, u, bc)
            hk.assign(h)
            uk.assign(u)

        collect_data()

        hn.assign(h)
        un.assign(u)


if __name__ == '__main__':
    shallow_water_1d(K=1, s=1, tau=0.01, mesh_size=100, T=1)
