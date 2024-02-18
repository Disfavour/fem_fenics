from fenics import *
import numpy as np
from scipy.constants import g
import matplotlib.pyplot as plt
from os.path import dirname, join

dpi=300

def plot_0(t, f, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, f)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$F$')
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_1(t, f, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, f)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$F1$')
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_2(t, f, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, f)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$F2$')
    plt.xlim(t[0], t[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()


set_log_level(LogLevel.WARNING)

def sw_1d(mesh_size=200, tau=0.005, s=1, vtkfname=None):
    vtkfile = File(vtkfname) if vtkfname is not None else None
    degree = 1
    T = 0.5

    m, E, me, Ee = [], [], [], []

    F_0, F_1, F_2 = [], [], []

    hl, hr = 10, 1
    h1 = 3.9618
    u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
    D1 = - (g * hl) ** 0.5
    D2 = u1 - (g*h1) ** 0.5
    D3 = u1 * h1 / (h1 - hr)
    domain_size = 5

    t = 0
    time_steps = np.linspace(t, T, round(T / tau) + 1)

    mesh = IntervalMesh(mesh_size, -domain_size, domain_size)

    H = FiniteElement("P", mesh.ufl_cell(), degree)
    U = FiniteElement("P", mesh.ufl_cell(), degree)
    W = FunctionSpace(mesh, MixedElement([H, U]))

    w, wn = Function(W), Function(W)
    h, u = split(w)
    hn, un = split(wn)


    w_exact = Expression(('x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))',
                          'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'),
                        g=g, hl=hl, hr=hr, h1=h1, u1=u1, D1=D1, D2=D2, D3=D3, t=t, degree=degree)

    hs = s*h + (1-s)*hn
    us = s*u + (1-s)*un

    F11 = (h-hn)/tau * dx
    F12 = (hs*us).dx(0) * dx

    F21 = (h*u-hn*un)/tau * dx
    F22 = (hs*us*us).dx(0) * dx
    F23 = g/2*(hs*hs).dx(0) * dx

    F1 = ((h-hn)/tau + (hs*us).dx(0)) * dx
    F2 = ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * dx
    F = F1 + F2


    # F = ((h-hn)/tau + (hs*us).dx(0)) * dx \
    #     + ((h*u-hn*un)/tau + (hs*us*us).dx(0) + g/2*(hs*hs).dx(0)) * dx
    
    m_eq = h * dx
    E_eq = 0.5*(h*u*u + g*h*h) * dx

    def collect_data():
        m.append(assemble(m_eq))
        E.append(assemble(E_eq))
        print(f'Time {t:>7.5f} m {m[-1]:>7.5f} E {E[-1]:>7.5f}')
        # print(assemble(F11), assemble(F12), assemble(F1))
        # print(assemble(F21), assemble(F22), assemble(F23), assemble(F2))
        # print(assemble(F))

        F_0.append(assemble(F))
        F_1.append(assemble(F1))
        F_2.append(assemble(F2))

        print(F_0[-1], F_1[-1], F_2[-1])
        print(D1*t, D2*t, D3*t)

        if vtkfile is not None:
            vtkfile << (w, t)

    # t = 0
    w.assign(project(w_exact, W))
    #collect_data()
    wn.assign(w)

    for t in time_steps[1:]:
        w_exact.t = t
        w.assign(project(w_exact, W))
        collect_data()
        wn.assign(w)

    return map(np.array, (time_steps, m, E, me, Ee, F_0, F_1, F_2))


if __name__ == '__main__':
    ms=1600
    tau=0.00125
    s=1

    base_dir = dirname(__file__)
    images_dir = join(base_dir, 'images')
    vtk_dir = join(base_dir, 'paraview')

    ts, m, E, me, Ee, F_0, F_1, F_2 = sw_1d(ms, tau, s)
    ts = ts[1:]

    

    plot_0(ts, F_0, join(images_dir, f'ms{ms}_tau{tau}_s{s}_F.png'))
    plot_1(ts, F_1, join(images_dir, f'ms{ms}_tau{tau}_s{s}_F1.png'))
    plot_2(ts, F_2, join(images_dir, f'ms{ms}_tau{tau}_s{s}_F2.png'))
