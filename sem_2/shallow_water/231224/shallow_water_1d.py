from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)

def sw_1d(mesh_size=200, tau=0.005, s=1):
    degree = 1
    T = 0.5

    m, E, me, Ee = [], [], [], []

    components = [[] for i in range(5)]
    components_exact = [[] for i in range(5)]

    hl, hr = 10, 1
    domain_size = 5
    h1 = 3.9618
    u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    D1 = - (g * hl) ** 0.5
    D2 = u1 - (g*h1) ** 0.5
    D3 = u1 * h1 / (h1 - hr)

    t = 0
    time_steps = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    x = mesh.coordinates().flatten()

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    bc = DirichletBC(W.sub(1), Constant(0), 'on_boundary')

    w, wn = Function(W), Function(W)
    ht, ut = TestFunctions(W)
    h, u = split(w)
    hn, un = split(wn)

    we, wne = Function(W), Function(W)
    he, ue = split(we)
    hne, une = split(wne)

    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    hse = s*he + (1-s)*hne
    use = s*ue + (1-s)*une

    # F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * ut*dx
    
    F = ((h-hn)/tau + (hs*us).dx(0)) * ht*dx \
        + ((h*u-hn*un)/tau + (hs*us*us).dx(0)) * ut*dx - g/2*(hs*hs) * ut.dx(0) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    m_eq_exact = he * dx
    E_eq_exact = 0.5*(he*ue*ue + g*he*he) * dx

    def collect_data():
        #m, E, me, Ee = map(assemble, (m_eq, E_eq, m_eq_exact, E_eq_exact))
        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        me.append(assemble(m_eq_exact))
        Ee.append(assemble(E_eq_exact))
        print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')

        components[0].append(assemble((h - hn)/tau * hs * dx))
        components[1].append(assemble((hs * us).dx(0) * hs * dx))
        components[2].append(assemble((h*u - hn*un)/tau * us * dx))
        components[3].append(assemble((hs*us*us).dx(0) * us * dx))
        components[4].append(assemble(g*hs*hs.dx(0) * us * dx))

        components_exact[0].append(assemble((he - hne)/tau * hse * dx))
        components_exact[1].append(assemble((hse * use).dx(0) * hse * dx))
        components_exact[2].append(assemble((he*ue - hne*une)/tau * use * dx))
        components_exact[3].append(assemble((hse*use*use).dx(0) * use * dx))
        components_exact[4].append(assemble(g*hse*hse.dx(0) * use * dx))
        

    # t = 0
    w.assign(project(Expression(('x[0] <= 0 ? hl : hr', '0'), hl=hl, hr=hr, degree=degree), W))

    w_exact.t = t
    we.assign(project(w_exact, W))

    collect_data()

    wn.assign(w)
    wne.assign(we)

    for t in time_steps[1:]:
        solve(F == 0, w, bc)

        w_exact.t = t
        we.assign(project(w_exact, W))

        collect_data()

        wn.assign(w)
        wne.assign(we)

    return map(np.array, (time_steps, m, E, me, Ee, components, components_exact))


if __name__ == '__main__':
    sw_1d()
    pass
