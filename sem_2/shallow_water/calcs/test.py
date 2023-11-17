from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join
from dolfin import *


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

set_log_level(LogLevel.WARNING)
#set_log_level(LogLevel.PROGRESS)


def shallow_water_1d(s=1, tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_size = 5
    h1 = 3.9618
    u1 = g * (h1*h1 - hr*hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    D1 = - (g*hl)**0.5
    D2 = u1 - (g*h1)**0.5
    D3 = u1*h1 / (h1 - hr)

    vtkfile = File(join(paraview, 'test.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bcs = [DirichletBC(W.sub(1), Constant(0), 'on_boundary')]

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

    F = ((h-hn)/tau + (hs()*us()).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs()*us()*us()).dx(0) + g/2*(hs()*hs()).dx(0)) * ut*dx
    
    # F = ((h-hn)/tau + (hs()*us()).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs()*us()*us()).dx(0) + g*hs()*(hs()).dx(0)) * ut*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

        vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0] < 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in ts[1:]:
        # J = derivative(F, w)

        # problem = NonlinearVariationalProblem(F, w , bcs, J)
        # solver = NonlinearVariationalSolver(problem)

        # prm = solver.parameters
        # prm['nonlinear_solver'] = 'snes'
        # prm['snes_solver']['line_search'] = 'bt'
        # prm['snes_solver']['linear_solver'] = 'mumps'
        # #prm[“snes_solver”][“linear_solver”] = “gmres”
        # #prm[“snes_solver”][“preconditioner”] = “hypre_amg”
        # #prm[“snes_solver”][“krylov_solver”][“nonzero_initial_guess”] = False
        # prm['snes_solver']['maximum_iterations'] = 1000
        # prm['snes_solver']['report'] = True

        # solver.solve()


        solve(F == 0, w, bcs, solver_parameters={'newton_solver':
                                                 {'absolute_tolerance': 1e-10,
                                                  'relative_tolerance': 1e-10,
                                                  'maximum_iterations': 10000,
                                                  'relaxation_parameter': 0.05}}) # 0.25

        collect_data()

        wn.assign(w)


if __name__ == '__main__':
    shallow_water_1d(s=0.5, tau=0.01, T=1, mesh_size=100)
