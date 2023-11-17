from fenics import *
import numpy as np

set_log_level(LogLevel.WARNING)


def sw_r_ru(ms=100, tau=0.01, degree_r=1, degree_w=1, sigma=1):
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
    W = VectorElement("P", mesh.ufl_cell(), degree_w)
    Q = FunctionSpace(mesh, MixedElement([R, W]))

    bcs = [
        DirichletBC(Q.sub(1).sub(0), Constant(0), CompiledSubDomain('on_boundary && (near(x[0], -size) || near(x[0], size))', size=size)),
        DirichletBC(Q.sub(1).sub(1), Constant(0), CompiledSubDomain('on_boundary && (near(x[1], -size) || near(x[1], size))', size=size))
    ]

    q, qn = Function(Q), Function(Q)
    rt, wt = TestFunctions(Q)
    r, w = split(q)
    rn, wn = split(qn)
    
    def rs():
        return sigma*r + (1-sigma)*rn
    
    def ws():
        return sigma*w + (1-sigma)*wn

    F = (r - rn)/tau*rt*dx + div(ws())*rt*dx \
        + dot((w - wn)/tau, wt)*dx + dot(div(rs()*outer(ws()/rs(), ws()/rs())), wt)*dx + dot(grad(a*rs()*rs()), wt)*dx

    m_eq = r * dx
    E_eq = (0.5*dot(w, w)/r + a*r*r) * dx

    def collect_data():
            m, E = map(assemble, (m_eq, E_eq))

            Es.append(E)

            rv = q.sub(0).compute_vertex_values()

            r_min.append(rv.min())
            r_max.append(rv.max())
            r_00.append(rv[idxs[0]])

            if np.isclose(t, moments).any():
                r_.append(rv[idxs])

            print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

    q.assign(project(Expression(('1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))', '0', '0'), alf=alf, bet=bet, degree=max(degree_r, degree_w)), Q))

    collect_data()

    qn.assign(q)

    for t in ts[1:]:
            solve(F == 0, q, bcs)

            collect_data()

            qn.assign(q)

    return map(np.array, (ts, Es, r_min, r_max, r_00, xs, r_))


if __name__ == '__main__':
    #sw_r_ru()
    pass
