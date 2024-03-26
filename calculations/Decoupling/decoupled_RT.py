from fenics import *
import numpy as np


def decoupled_RT(K=1, degree=2, mesh_size=100, tau=0.005, T=5, t_fig4=(1, 2, 3, 4, 5)):
    set_log_level(LogLevel.WARNING)

    ts = []
    Es = []
    xs = np.linspace(0, 5, int(T / tau) + 1)
    rs = []

    t = 0
    time_steps = [(i + 1) * tau for i in range(int(T / tau))]

    domain_size = 5

    alf = 2
    bet = 20

    tau_ = Constant(tau)
    a = Constant(1)
    gam = Constant(1.4)

    mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size,
                         "right")

    V = FunctionSpace(mesh, "RT", 1)
    S = FunctionSpace(mesh, "P", degree)

    bc = DirichletBC(V, Constant((0, 0)), 'on_boundary')

    ut, rt = TestFunction(V), TestFunction(S)
    u_, r_ = TrialFunction(V), TrialFunction(S)
    u, r = Function(V), Function(S)
    un, rn = Function(V), Function(S)
    uk, rk = Function(V), Function(S)

    F1 = (r_ - rn) / tau_ * rt * dx \
         - dot(r_ * uk, grad(rt)) * dx

    a1, L1 = lhs(F1), rhs(F1)

    F2 = dot((r * u_ - rn * un) / tau_, ut) * dx \
         - inner(outer(r * uk, u_), nabla_grad(ut)) * dx \
         + dot(grad(a * r ** gam), ut) * dx

    a2, L2 = lhs(F2), rhs(F2)

    m_eq = r * dx
    E_eq = (r * dot(u, u) / 2 + a * r ** gam / (gam - 1)) * dx

    # t = 0
    u.assign(project(Constant((0, 0)), V))
    r.assign(project(Expression("1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=degree), S))

    m, E = assemble(m_eq), assemble(E_eq)
    print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

    ts.append(t)
    Es.append(E)

    un.assign(u)
    rn.assign(r)

    for t in time_steps:
        uk.assign(un)
        rk.assign(rn)

        for k in range(K):
            solve(a1 == L1, r)
            solve(a2 == L2, u, bc)
            uk.assign(u)
            rk.assign(r)

        m, E = assemble(m_eq), assemble(E_eq)
        print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

        ts.append(t)
        Es.append(E)

        for t4 in t_fig4:
            if abs(t4 - t) < tau / 2:
                ros_cur = []
                for x1 in xs:
                    ros_cur.append(r((x1, 0)))
                rs.append(np.array(ros_cur))

        un.assign(u)
        rn.assign(r)

    ts = np.array(ts)
    Es = np.array(Es)

    np.save(f"data/dec_RT_K{K}_d{degree}_ms{mesh_size}_tau{tau}.npy", np.array((ts, Es, xs, *rs)).T)


if __name__ == "__main__":
    pass
