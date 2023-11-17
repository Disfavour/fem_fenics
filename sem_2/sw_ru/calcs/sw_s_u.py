from fenics import *
import numpy as np

set_log_level(LogLevel.WARNING)


def sw_s_u(ms=100, tau=0.01, degree_s=1, degree_u=1, sigma=1):
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

    S = FiniteElement('P', mesh.ufl_cell(), degree_s)
    U = VectorElement("P", mesh.ufl_cell(), degree_u)
    Q = FunctionSpace(mesh, MixedElement([S, U]))

    bcs = [
        DirichletBC(Q.sub(1).sub(0), Constant(0), CompiledSubDomain('on_boundary && (near(x[0], -size) || near(x[0], size))', size=size)),
        DirichletBC(Q.sub(1).sub(1), Constant(0), CompiledSubDomain('on_boundary && (near(x[1], -size) || near(x[1], size))', size=size))
    ]

    q, qn = Function(Q), Function(Q)
    st, ut = TestFunctions(Q)
    s, u = split(q)
    sn, un = split(qn)
    
    def ss():
        return sigma*s + (1-sigma)*sn
    
    def us():
        return sigma*u + (1-sigma)*un

    F = (s*s - sn*sn)/tau*st*dx + div(ss()*ss()*us())*st*dx \
        + dot((s*s*u - sn*sn*un)/tau, ut)*dx + dot(div(ss()*ss()*outer(us(), us())), ut)*dx + dot(grad(a*ss() ** 4), ut)*dx

    m_eq = s*s * dx
    E_eq = (0.5*s*s*dot(u, u) + a*s ** 4) * dx

    def collect_data():
            m, E = map(assemble, (m_eq, E_eq))

            Es.append(E)

            rv = q.sub(0).compute_vertex_values() ** 2

            r_min.append(rv.min())
            r_max.append(rv.max())
            r_00.append(rv[idxs[0]])

            if np.isclose(t, moments).any():
                r_.append(rv[idxs])

            print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

    q.assign(project(Expression(('sqrt(1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])))', '0', '0'), alf=alf, bet=bet, degree=max(degree_s, degree_u)), Q))

    collect_data()

    qn.assign(q)

    for t in ts[1:]:
            solve(F == 0, q, bcs)

            collect_data()

            qn.assign(q)

    return map(np.array, (ts, Es, r_min, r_max, r_00, xs, r_))


if __name__ == '__main__':
    sw_s_u(ms=10, tau=0.1, sigma=0.5)
    pass
