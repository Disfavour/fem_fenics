from fenics import *
import numpy as np
from os.path import dirname, join

set_log_level(LogLevel.WARNING)

base_dir = dirname(dirname(__file__))
data = join(base_dir, 'data')
paraview = join(base_dir, 'paraview')


def symmetric(tau=0.01, mesh_size=100, vtkfile=None, verbose=False):
    fname = f'sym_tau{tau}_ms{mesh_size}'

    ts = []
    ms = []
    I0s, I1s = [], []
    dI0s, dI1s = [], []
    Es = []
    rs = []
    r_00 = []
    r_min = []
    r_max = []

    degree = 1
    t, T = 0, 5
    time_steps = [(i + 1) * tau for i in range(int(T / tau))]

    t_fig4 = (1, 2, 3, 4, 5)
    xs = np.linspace(0, 5, len(time_steps) + 1)

    domain_size = 5

    alf = 2
    bet = 20

    tau_ = Constant(tau)
    a = Constant(1)
    gam = Constant(2)

    mesh = RectangleMesh(Point(-domain_size, -domain_size), Point(domain_size, domain_size), mesh_size, mesh_size,
                         "right")
    n = FacetNormal(mesh)

    W = VectorElement("P", mesh.ufl_cell(), degree)
    S = FiniteElement("P", mesh.ufl_cell(), degree)
    V = FunctionSpace(mesh, MixedElement([W, S]))

    R = FunctionSpace(mesh, "P", degree)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_x = CompiledSubDomain("on_boundary && (near(x[0], -size) || near(x[0], size))", size=domain_size)
    subdomain_y = CompiledSubDomain("on_boundary && (near(x[1], -size) || near(x[1], size))", size=domain_size)
    subdomain_x.mark(boundaries, 1)
    subdomain_y.mark(boundaries, 2)

    bcs = [DirichletBC(V.sub(0).sub(0), Constant(0), boundaries, 1),
           DirichletBC(V.sub(0).sub(1), Constant(0), boundaries, 2)]

    wt, st = TestFunctions(V)
    v = Function(V)
    vn = Function(V)
    w, s = split(v)
    wn, sn = split(vn)

    def wm():
        return (w + wn) / 2

    def sm():
        return (s + sn) / 2

    def um():
        return wm() / sm()

    F = (s - sn) / tau_ * st * dx \
        + (div(um() * sm()) + dot(um(), grad(sm()))) / 2 * st * dx \
        + dot((w - wn) / tau_, wt) * dx \
        + dot(div(outer(um(), wm())) + dot(um(), grad(wm())) / 2, wt) * dx \
        + dot(1 / sm() * grad(a * (sm() ** 2) ** gam), wt) * dx

    m_eq = s ** 2 * dx
    E_eq = (dot(w, w) / 2 + a * (s ** 2) ** gam) / (gam - 1) * dx
    I_eq0 = s * w[0] * dx
    I_eq1 = s * w[1] * dx
    dI_eq0 = tau * a * (sm() ** 2) ** gam * n[0] * ds
    dI_eq1 = tau * a * (sm() ** 2) ** gam * n[1] * ds

    # t = 0
    v.assign(project(
        Expression(("0", "0", "pow(1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1])), 0.5)"), alf=alf, bet=bet, degree=degree),
        V))

    m, E = assemble(m_eq), assemble(E_eq)
    I0, I1 = assemble(I_eq0), assemble(I_eq1)
    dI0, dI1 = assemble(dI_eq0), assemble(dI_eq1)

    ts.append(t)
    ms.append(m)
    Es.append(E)
    I0s.append(I0)
    I1s.append(I1)
    dI0s.append(dI0)
    dI1s.append(dI1)

    w_, s_ = v.split()
    r_ = project(s_ ** 2, R)

    r_00.append(r_(0, 0))
    r_min.append(r_.compute_vertex_values().min())
    r_max.append(r_.compute_vertex_values().max())

    if verbose:
        print(f"conservation: t {t:>5.3f}\t M {m:>15.10f}\t E {E:>15.10f}\t"
              f"I ({I0:>15.10f}, {I1:>15.10f})\t I_0 - sum(dI) ({I0s[0] - sum(dI0s):>15.10f}, {I1s[0] - sum(dI1s):>15.10f})")

    if vtkfile is not None:
        r_.rename('r', 'label')
        vtkfile << (r_, t)

    vn.assign(v)

    for t in time_steps:
        solve(F == 0, v, bcs)

        m, E = assemble(m_eq), assemble(E_eq)
        I0, I1 = assemble(I_eq0), assemble(I_eq1)
        dI0, dI1 = assemble(dI_eq0), assemble(dI_eq1)

        ts.append(t)
        ms.append(m)
        Es.append(E)
        I0s.append(I0)
        I1s.append(I1)
        dI0s.append(dI0)
        dI1s.append(dI1)

        w_, s_ = v.split()
        r_ = project(s_ ** 2, R)

        r_00.append(r_(0, 0))
        r_min.append(r_.compute_vertex_values().min())
        r_max.append(r_.compute_vertex_values().max())

        if verbose:
            print(f"conservation: t {t:>5.3f}\t M {m:>15.10f}\t E {E:>15.10f}\t"
                  f"I ({I0:>15.10f}, {I1:>15.10f})\t I_0 - sum(dI) ({I0s[0] - sum(dI0s):>15.10f}, {I1s[0] - sum(dI1s):>15.10f})")

        for t4 in t_fig4:
            if abs(t4 - t) < tau / 2:
                ros_cur = []
                for x1 in xs:
                    ros_cur.append(r_((x1, 0)))
                rs.append(np.array(ros_cur))

        if vtkfile is not None:
            r_.rename('r', 'label')
            vtkfile << (r_, t)

        vn.assign(v)

    np.save(join(data, fname + '.npy'), np.array((ts, ms, Es, I0s, I1s, dI0s, dI1s, r_00, r_min, r_max, xs, *rs)))


if __name__ == '__main__':
    symmetric(verbose=True)
    #symmetric(tau=0.01, vtkfile=File(join(paraview, 'sym_t1' + '.pvd')), verbose=True)
