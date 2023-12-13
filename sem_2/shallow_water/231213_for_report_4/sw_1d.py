from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)

def sw_1d(s, tau, mesh_size, T, dif_t, critical_time, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None
    h_numerical_list = []
    u_numerical_list = []
    h_exact_list = []
    u_exact_list = []

    ms = []
    Es = []

    hl, hr = 10, 1
    domain_size = 5
    h1 = 3.9618
    u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    D1 = - (g * hl) ** 0.5
    D2 = u1 - (g*h1) ** 0.5
    D3 = u1 * h1 / (h1 - hr)

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

    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)

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

        if t <= critical_time:
            w_exact.t = t
            we.assign(project(w_exact, W))

            if np.isclose(t, dif_t).any():
                h_numerical_list.append(w.sub(0).compute_vertex_values())
                u_numerical_list.append(w.sub(1).compute_vertex_values())
                h_exact_list.append(we.sub(0).compute_vertex_values())
                u_exact_list.append(we.sub(1).compute_vertex_values())

        else:
            if np.isclose(t, dif_t).any():
                h_numerical_list.append(w.sub(0).compute_vertex_values())
                u_numerical_list.append(w.sub(1).compute_vertex_values())
                h_exact_list.append(np.full(x.shape, np.nan))
                u_exact_list.append(np.full(x.shape, np.nan))
        
        if vtkfile is not None:
            vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0] < 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in time_steps[1:]:
        solve(F == 0, w, bc)

        collect_data()

        wn.assign(w)

    return map(np.array, (x, h_numerical_list, u_numerical_list, h_exact_list, u_exact_list, time_steps, Es))


if __name__ == '__main__':
    # shallow_water_1d(1, 0.01, 100, 1, [], 0)
    # sw_1d(s, tau, mesh_size, T, dif_t, critical_time, degree=1, vtkfname=None)
    pass
