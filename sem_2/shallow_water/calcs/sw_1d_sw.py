from fenics import *
import numpy as np
from scipy.constants import g


set_log_level(LogLevel.WARNING)


def shallow_water_1d(sigma=1, tau=0.01, mesh_size=100, T=1, degree=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None

    hl, hr = 10, 1
    domain_size = 5

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    S = FiniteElement('P', mesh.ufl_cell(), degree)
    W = FiniteElement('P', mesh.ufl_cell(), degree)
    Q = FunctionSpace(mesh, MixedElement([S, W]))

    bc = DirichletBC(Q.sub(1), Constant(0), 'on_boundary')

    q, qn = Function(Q), Function(Q)
    st, wt = TestFunctions(Q)
    s, w = split(q)
    sn, wn = split(qn)

    q.rename('(s, w)', 'shallow water 1d')
    
    ss = sigma*s + (1-sigma)*sn
    ws = sigma*w + (1-sigma)*wn
    
    F = (s - sn)/tau*st*dx \
        + 0.5/ss*(ss*ws).dx(0) * st * dx \
        + (s*w - sn*wn)/tau * wt*dx \
        + (ws*ws).dx(0)*wt*dx \
        + 2*g*ss**3*ss.dx(0) *wt * dx
    
    # F = (s - sn)/tau * st * dx \
    #     + (ws.dx(0) + ws/ss*ss.dx(0)) / 2 * st * dx \
    #     + (w - wn)/tau * wt * dx \
    #     + 0.5*((ws/ss*ws).dx(0) + ws/ss*ws.dx(0)) * wt * dx \
    #     + 2*g*ss**2*ss.dx(0) * wt * dx
    
    m_eq = s*s * dx
    E_eq = 0.5*(w*w + g*s**4) * dx
    
    def collect_data():
        m, E = map(assemble, (m_eq, E_eq))

        if vtkfile is not None:
            vtkfile << (q, t)

        print(f'Time {t:>7.5f} m {m:>7.8f} E {E:>7.8f}')

    # t = 0
    q.assign(project(Expression(('x[0] < 0 ? sqrt(hl) : sqrt(hr)', '0'), hl=hl, hr=hr, degree=degree), Q))

    collect_data()

    qn.assign(q)

    for t in ts[1:]:
        solve(F == 0, q, bc)

        collect_data()

        qn.assign(q)


if __name__ == '__main__':
    from os.path import dirname, join
    base_dir = dirname(dirname(__file__))
    paraview = join(base_dir, 'paraview')
    vtkfname = join(paraview, 'numerical_sw.pvd')
    shallow_water_1d(sigma=0.5, tau=0.005, mesh_size=200, vtkfname=None)
