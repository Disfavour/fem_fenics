from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)

def sw_2d_v1(s, tau, mesh_size, T, dif_t, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None

    ms = []
    Es = []
    hss = []
    hps = []

    hl, hr = 3, 1
    domain_size_x = 5
    domain_size_y = 1
    mesh_size_y = mesh_size * domain_size_y // domain_size_x

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = RectangleMesh(Point(-domain_size_x, -domain_size_y), Point(domain_size_x, domain_size_y), mesh_size, mesh_size_y)

    x = mesh.coordinates()
    idxs = np.where(x[:, 1] == 0)
    x = x[idxs][:, 0]
    #print(x)

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = VectorElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bcs = [
        DirichletBC(W.sub(1).sub(0), Constant(0), CompiledSubDomain("on_boundary && (near(x[0], -xsize) || near(x[0], xsize))", xsize=domain_size_x)),
        DirichletBC(W.sub(1).sub(1), Constant(0), CompiledSubDomain("on_boundary && (near(x[1], -ysize) || near(x[1], ysize))", ysize=domain_size_y)),
    ]

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    F = ((h-hn)/tau + div(hs*us)) * ht * dx \
        + dot(((h*u-hn*un)/tau + div(hs*outer(us, us)) + g*hs*grad(hs)), ut) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*dot(u, u) + g*h*h) * dx

    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))
        print(f'Time {t:>7.5f} m {m:>7.5f} E {E:>7.5f}')

        ms.append(m)
        Es.append(E)

        if np.isclose(t, dif_t).any():
            hss.append(w.sub(0).compute_vertex_values()[idxs])
        
        if vtkfile is not None:
            vtkfile << (w, t)

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0', '0'), hl=hl, hr=hr, degree=degree), W))

    collect_data()

    wn.assign(w)

    for t in ts[1:]:
        solve(F == 0, w, bcs)
        collect_data()
        wn.assign(w)

    return map(np.array, (ts, ms, Es, x, hss, hps))


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(__file__)
    paraview = join(base_dir, 'paraview')
    vtkfname = join(paraview, 'sw_2d_hl3_s1_tau0.005_m200.pvd')
    #sw_2d(0.5, 0.005, 200, 2, degree=1, vtkfname=None)
    pass
