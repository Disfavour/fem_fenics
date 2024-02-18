from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)

def sw_1d(hl, s, tau, mesh_size, T, dif_t, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None
    h_numerical_list = []
    u_numerical_list = []
    h_exact_list = []
    u_exact_list = []

    ms = []
    Es = []

    hr = 1
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
    we = Function(W)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * ut*dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))
        print(f'Time {t:>7.5f} m {m:>7.5f} E {E:>7.5f}')

        ms.append(m)
        Es.append(E)

        if np.isclose(t, dif_t).any():
            h_numerical_list.append(w.sub(0).compute_vertex_values())
            u_numerical_list.append(w.sub(1).compute_vertex_values())
        
        if vtkfile is not None:
            vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in time_steps[1:]:
        solve(F == 0, w, bc)

        collect_data()

        wn.assign(w)

    return map(np.array, (x, h_numerical_list, u_numerical_list, time_steps, Es))


if __name__ == '__main__':
    # shallow_water_1d(1, 0.01, 100, 1, [], 0)
    # sw_1d(s, tau, mesh_size, T, dif_t, critical_time, degree=1, vtkfname=None)
    pass
