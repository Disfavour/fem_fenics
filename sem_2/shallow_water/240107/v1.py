from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)

def sw_1d(mesh_size=200, tau=0.005, s=1):
    degree = 1
    T = 2

    m, E = [], []

    hl, hr = 3, 1
    domain_size = 5

    t = 0
    time_steps = np.linspace(t, T, round(T / tau) + 1)

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
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g*hs*hs.dx(0)) * ut*dx
    
    # F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0)) * ut*dx - g/2*(hs*hs) * ut.dx(0) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    def collect_data():
        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))
    collect_data()
    wn.assign(w)

    for t in time_steps[1:]:
        solve(F == 0, w, bc)
        collect_data()
        wn.assign(w)

    return map(np.array, (time_steps, m, E))


if __name__ == '__main__':
    sw_1d(s=0.5)
    pass
