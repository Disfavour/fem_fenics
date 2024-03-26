from fenics import *
import numpy as np


def seidel(K=1, degree=1, mesh_size=100, tau=0.005, T=5, t_fig4=(1, 2, 3, 4, 5)):
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

    U = FunctionSpace(mesh, "P", degree)
    V = FunctionSpace(mesh, "P", degree)
    R = FunctionSpace(mesh, "P", degree)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
    subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
    subdomain_x.mark(boundaries, 1)
    subdomain_y.mark(boundaries, 2)

    bc_u = DirichletBC(U, Constant(0), boundaries, 1)
    bc_v = DirichletBC(V, Constant(0), boundaries, 2)

    ut, vt, rt = TestFunction(U), TestFunction(V), TestFunction(R)
    u_, v_, r_ = TrialFunction(U), TrialFunction(V), TrialFunction(R)
    u, v, r = Function(U), Function(V), Function(R)
    un, vn, rn = Function(U), Function(V), Function(R)
    uk, vk, rk = Function(U), Function(V), Function(R)

    F1 = (r_ - rn) / tau_ * rt * dx \
         + (r_ * uk).dx(0) * rt * dx + (r_ * vk).dx(1) * rt * dx

    a1, L1 = lhs(F1), rhs(F1)

    F2 = (r * u_ - rn * un) / tau_ * ut * dx \
         + (r * u_ * uk).dx(0) * ut * dx + (r * u_ * vk).dx(1) * ut * dx \
         + a * gam * r ** (gam - 1) * r.dx(0) * ut * dx

    a2, L2 = lhs(F2), rhs(F2)

    F3 = (r * v_ - rn * vn) / tau_ * vt * dx \
         + (r * v_ * u).dx(0) * vt * dx + (r * v_ * vk).dx(1) * vt * dx \
         + a * gam * r ** (gam - 1) * r.dx(1) * vt * dx

    a3, L3 = lhs(F3), rhs(F3)

    m_eq = r * dx
    E_eq = 0.5 * r * (u * u + v * v) * dx + 1. / (gam - 1) * r ** gam * dx

    # t = 0
    u.assign(project(Constant(0), U))
    v.assign(project(Constant(0), V))
    r.assign(project(Expression("1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))", alf=alf, bet=bet, degree=degree), R))

    m, E = assemble(m_eq), assemble(E_eq)
    print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

    ts.append(t)
    Es.append(E)

    un.assign(u)
    vn.assign(v)
    rn.assign(r)

    for t in time_steps:
        uk.assign(un)
        vk.assign(vn)
        rk.assign(rn)

        for k in range(K):
            solve(a1 == L1, r)
            solve(a2 == L2, u, bc_u)
            solve(a3 == L3, v, bc_v)

            uk.assign(u)
            vk.assign(v)
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
        vn.assign(v)
        rn.assign(r)

    ts = np.array(ts)
    Es = np.array(Es)

    np.save(f"data/sei_K{K}_d{degree}_ms{mesh_size}_tau{tau}.npy", np.array((ts, Es, xs, *rs)).T)


if __name__ == "__main__":
    seidel(tau=0.02, T=1)
    pass
