from fenics import *
import numpy as np


def nonlinear(degree=1, mesh_size=100, tau=0.005, T=5, t_fig4=(1, 2, 3, 4, 5)):
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

    mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size, "right")

    V = VectorElement("P", mesh.ufl_cell(), degree)
    S = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([V, S]))

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
    subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
    subdomain_x.mark(boundaries, 1)
    subdomain_y.mark(boundaries, 2)

    bcs = [DirichletBC(W.sub(0).sub(0), Constant(0), boundaries, 1),
           DirichletBC(W.sub(0).sub(1), Constant(0), boundaries, 2)]

    ut, rt = TestFunctions(W)
    w = Function(W)
    wn = Function(W)
    u, r = split(w)
    un, rn = split(wn)

    F = (r-rn)/tau_ * rt * dx\
        - dot(r*u, grad(rt)) * dx\
        + dot((r*u-rn*un)/tau_, ut) * dx\
        - inner(outer(r*u, u), nabla_grad(ut)) * dx\
        + dot(grad(a * r**gam), ut) * dx

    m_eq = r * dx
    E_eq = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx

    # t = 0
    w.assign(project(Expression(("0", "0", "1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))"), alf=alf, bet=bet, degree=degree), W))

    m, E = assemble(m_eq), assemble(E_eq)
    print(f"conservation: t {t:.3f}\t M {m:.10f}\t E {E:.10f}")

    ts.append(t)
    Es.append(E)

    wn.assign(w)

    for t in time_steps:
        solve(F == 0, w, bcs)

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

        wn.assign(w)

    ts = np.array(ts)
    Es = np.array(Es)

    np.save(f"data/nonlinear_d{degree}_ms{mesh_size}_tau{tau}.npy", np.array((ts, Es, xs, *rs)).T)


if __name__ == "__main__":
    pass
