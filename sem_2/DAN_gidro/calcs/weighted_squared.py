from fenics import *
import numpy as np
from os.path import dirname, join

set_log_level(LogLevel.WARNING)

base_dir = dirname(dirname(__file__))
data = join(base_dir, 'data')
paraview = join(base_dir, 'paraview')


def weighted(s=1, tau=0.01, mesh_size=100, vtkfile=None, verbose=False):
    fname = f'w_s{s}_tau{tau}_ms{mesh_size}'

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

    def rm():
        return (s*r + (1-s)*rn)

    def um():
        return (s*u + (1-s)*un)

    F = (r-rn)/tau_ * rt * dx \
        - dot(rm()*um(), grad(rt)) * dx \
        + dot((r*u-rn*un)/tau_, ut) * dx \
        - inner(outer(rm()*um(), um()), nabla_grad(ut)) * dx \
        + dot(grad(a * rm()**gam), ut) * dx

    m_eq = r * dx
    E_eq = (r*dot(u, u)/2 + a*r**gam/(gam-1)) * dx
    I_eq0 = r * u[0] * dx
    I_eq1 = r * u[1] * dx
    dI_eq0 = tau * a * r ** gam * n[0] * ds
    dI_eq1 = tau * a * r ** gam * n[1] * ds

    # t = 0
    w.assign(project(Expression(("0", "0", "1 + alf*exp(-bet*(x[0]*x[0]+x[1]*x[1]))"), alf=alf, bet=bet, degree=degree), W))

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

    u_, r_ = w.split()

    r_00.append(r_(0, 0))
    r_min.append(r_.compute_vertex_values().min())
    r_max.append(r_.compute_vertex_values().max())

    if verbose:
        print(f"conservation: t {t:>5.3f}\t M {m:>15.10f}\t E {E:>15.10f}\t"
              f"I ({I0:>15.10f}, {I1:>15.10f})\t I_0 - sum(dI) ({I0s[0] - sum(dI0s):>15.10f}, {I1s[0] - sum(dI1s):>15.10f})")

    if vtkfile is not None:
        r_.rename('r', 'label')
        vtkfile << (r_, t)

    wn.assign(w)

    for t in time_steps:
        solve(F == 0, w, bcs)

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

        u_, r_ = w.split()

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
                    ros_cur.append(r((x1, 0)))
                rs.append(np.array(ros_cur))

        if vtkfile is not None:
            r_.rename('r', 'label')
            vtkfile << (r_, t)

        wn.assign(w)

    np.save(join(data, fname + '.npy'), np.array((ts, ms, Es, I0s, I1s, dI0s, dI1s, r_00, r_min, r_max, xs, *rs)))


if __name__ == '__main__':
    weighted(0, tau=0.005, verbose=True)
    #weighted(tau=0.01, vtkfile=File(join(paraview, 'nsym_t1' + '.pvd')), verbose=True)
