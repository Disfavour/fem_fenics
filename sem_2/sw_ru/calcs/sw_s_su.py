from fenics import *
import numpy as np

set_log_level(LogLevel.WARNING)


def sw_s_su(ms=100, tau=0.01, degree_s=1, degree_w=1, sigma=1):
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
    W = VectorElement("P", mesh.ufl_cell(), degree_w)
    Q = FunctionSpace(mesh, MixedElement([S, W]))

    bcs = [
        DirichletBC(Q.sub(1).sub(0), Constant(0), CompiledSubDomain('on_boundary && (near(x[0], -size) || near(x[0], size))', size=size)),
        DirichletBC(Q.sub(1).sub(1), Constant(0), CompiledSubDomain('on_boundary && (near(x[1], -size) || near(x[1], size))', size=size))
    ]

    q, qn = Function(Q), Function(Q)
    st, wt = TestFunctions(Q)
    s, w = split(q)
    sn, wn = split(qn)
    
    def ss():
        return sigma*s + (1-sigma)*sn
    
    def ws():
        return sigma*w + (1-sigma)*wn

    # F = (s*s - sn*sn)/tau*st*dx + div(ss()*ws())*st*dx \
    #     + dot((s*w - sn*wn)/tau, wt)*dx + dot(div(ss()*ss()*outer(ws()/ss(), ws()/ss())), wt)*dx + dot(grad(a*ss() ** 4), wt)*dx
    
    # F = (s - sn)/tau*st*dx \
    #     + 0.5/ss()*div(ss()*ws())*st*dx \
    #     + dot((s*w - sn*wn)/tau, wt)*dx \
    #     + dot(div(outer(ws(), ws())), wt)*dx \
    #     + dot(a*4*ss()**3*grad(ss()), wt)*dx

    F = (s - sn) / tau * st * dx \
        + (div(ws()) + dot(ws()/ss(), grad(ss()))) / 2 * st * dx \
        + dot((w - wn) / tau, wt) * dx \
        + 0.5*dot(div(outer(ws()/ss(), ws())) + dot(ws()/ss(), grad(ws())), wt) * dx \
        + dot(1 / ss() * grad(a*ss() ** 4), wt) * dx
    
    # dot((w - wn) / tau, wt) * dx

    m_eq = s*s * dx
    E_eq = (0.5*dot(w, w) + a*s ** 4) * dx

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

    q.assign(project(Expression(('sqrt(1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])))', '0', '0'), alf=alf, bet=bet, degree=max(degree_s, degree_w)), Q))

    collect_data()

    qn.assign(q)

    for t in ts[1:]:
            solve(F == 0, q, bcs)

            collect_data()

            qn.assign(q)

    return map(np.array, (ts, Es, r_min, r_max, r_00, xs, r_))


if __name__ == '__main__':
    sw_s_su(sigma=0.5)
    pass
