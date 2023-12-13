from fenics import *
import numpy as np

set_log_level(LogLevel.WARNING)


def sw_r_u(ms=100, tau=0.01, degree_r=1, degree_u=1, sigma=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None

    size = 5
    alf = 2
    bet = 20

    t, T = 0, 5
    ts = np.linspace(t, T, round(T / tau) + 1)

    tau = Constant(tau)
    a = Constant(1)  # g/2

    moments = (1, 2, 3)
    Es = []
    r_min, r_max, r_00 = [], [], []
    r_ = []

    mesh = RectangleMesh(Point(-size, -size), Point(size, size), ms, ms)

    x = mesh.coordinates()
    idxs = np.where((x[:, 1] == 0) & (x[:, 0] >= 0))[0]
    xs = x[idxs][:, 0]

    R = FiniteElement('P', mesh.ufl_cell(), degree_r)
    U = VectorElement("P", mesh.ufl_cell(), degree_u)
    Q = FunctionSpace(mesh, MixedElement([R, U]))

    bcs = [
        DirichletBC(Q.sub(1).sub(0), Constant(0), CompiledSubDomain('on_boundary && (near(x[0], -size) || near(x[0], size))', size=size)),
        DirichletBC(Q.sub(1).sub(1), Constant(0), CompiledSubDomain('on_boundary && (near(x[1], -size) || near(x[1], size))', size=size))
    ]

    q, qn = Function(Q), Function(Q)
    rt, ut = TestFunctions(Q)
    r, u = split(q)
    rn, un = split(qn)

    rs = sigma*r + (1-sigma)*rn
    us = sigma*u + (1-sigma)*un

    F = (r - rn)/tau*rt*dx + div(rs*us)*rt*dx \
        + dot((r*u - rn*un)/tau, ut)*dx + dot(div(rs*outer(us, us)), ut)*dx + dot(grad(a*rs*rs), ut)*dx

    m_eq = r * dx
    E_eq = (0.5*r*dot(u, u) + a*r*r) * dx

    def collect_data():
            m, E = map(assemble, (m_eq, E_eq))

            Es.append(E)

            rv = q.sub(0).compute_vertex_values()

            r_min.append(rv.min())
            r_max.append(rv.max())
            r_00.append(rv[idxs[0]])

            if np.isclose(t, moments).any():
                r_.append(rv[idxs])

            if vtkfile is not None:
                vtkfile << (q, t)

            print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')
            print(r_min[-1])

    q.assign(project(Expression(('0.15 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))', '0', '0'), alf=alf, bet=bet, degree=max(degree_r, degree_u)), Q))

    collect_data()

    qn.assign(q)

    for t in ts[1:]:
            solve(F == 0, q, bcs)

            collect_data()

            qn.assign(q)

    return map(np.array, (ts, Es, r_min, r_max, r_00, xs, r_))


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(dirname(__file__))
    paraview = join(base_dir, 'paraview')
    vtkfname=join(paraview, 'dry_bottom_s0.5_r_u.pvd')
    sw_r_u(sigma=0.5, vtkfname=vtkfname)
    pass
