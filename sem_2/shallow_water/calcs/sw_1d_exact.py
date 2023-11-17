from fenics import *
import numpy as np
from scipy.constants import g
from os.path import dirname, join


base_dir = dirname(dirname(__file__))
paraview = join(base_dir, 'paraview')

set_log_level(LogLevel.WARNING)

def shallow_water_1d(tau=0.01, mesh_size=100, T=1, degree=1):
    hl, hr = 10, 1
    domain_size = 5
    h1 = 3.9618
    u1 = g * (h1*h1 - hr*hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    D1 = - (g*hl)**0.5
    D2 = u1 - (g*h1)**0.5
    D3 = u1*h1 / (h1 - hr)

    vtkfile = File(join(paraview, 'exact.pvd'))

    t = 0
    ts = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    w = Function(W)
    w.rename('(h, u)', 'shallow water 1d')

    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)
    
    
    for t in ts:
        print(f'Time {t}')

        w_exact.t = t
        w.assign(project(w_exact, W))
        vtkfile << (w, t)



if __name__ == '__main__':
    shallow_water_1d(mesh_size=1000)
