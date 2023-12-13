import matplotlib.pyplot as plt
import numpy as np
from os.path import join


l = ('-', '--', ':') # '-.'
c = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
lw = [0.9, 1.2, 1.5]
dpi = 600


def base_plot(x, f_numerical, f_exact):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(x, f_numerical, '-', label=r"numerical")
    plt.plot(x, f_exact, ':', label=r"exact")

    if not np.isnan(f_exact).all():
        plt.legend()

    plt.grid()
    plt.xlim(x[0], x[-1])
    plt.xlabel(r'$x$')


def h(x, h_numerical, h_exact, dir, filename):
    base_plot(x, h_numerical, h_exact)
    
    plt.ylabel(r'$h$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()


def u(x, u_numerical, u_exact, dir, filename):
    base_plot(x, u_numerical, u_exact)

    plt.ylabel(r'$u$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()


def at_different_moments(x, f_numerical, f_exact, time_moments):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for i, t in enumerate(time_moments):
        plt.plot(x, f_numerical[i], f'-{c[i]}', label=fr'$t={t}$')
        plt.plot(x, f_exact[i], f':{c[i]}')

    plt.legend()
    plt.grid()
    plt.xlim(x[0], x[-1])
    plt.xlabel(r'$x$')


def h_at_different_moments(x, h_numerical, h_exact, time_moments, dir, filename):
    at_different_moments(x, h_numerical, h_exact, time_moments)

    plt.ylabel(r'$h$')

    for i in ['pdf', 'png']:
        plt.savefig(join(dir, i, f'{filename}.{i}'))
    
    plt.close()


def u_at_different_moments(x, u_numerical, u_exact, time_moments, dir, filename):
    at_different_moments(x, u_numerical, u_exact, time_moments)

    plt.ylabel(r'$u$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()


def L2(t, f):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    
    plt.plot(t, f)

    plt.grid()
    plt.xlim(t[0], t[-1])
    plt.xlabel(r'$t$')


def h_L2(t, h, dir, filename):
    L2(t, h)

    plt.ylabel(r'$|| h_{exact} - h_{numerical} ||_2$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()


def u_L2(t, u, dir, filename):
    L2(t, u)

    plt.ylabel(r'$|| u_{exact} - u_{numerical} ||_2$')

    for i in ['pdf', 'png']:
        plt.savefig(join(dir, i, f'{filename}.{i}'))
    
    plt.close()


def E_different_s(t, Es, sigmas, dir, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    
    for i, (s, E) in enumerate(zip(sigmas, Es)):
        plt.plot(t, E, c[i], label=fr'$\sigma={s}$')

    plt.legend()
    plt.grid()
    plt.xlim(t[0], t[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()

def E_different_tau(ts, Es, taus, dir, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    
    for i, (tau, t, E) in enumerate(zip(taus, ts, Es)):
        plt.plot(t, E, c[i], label=fr'$\tau={tau}$')

    plt.legend()
    plt.grid()
    plt.xlim(t[0], t[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()

def delta_different_s(t, deltas, sigmas, dir, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    
    for i, (s, delta) in enumerate(zip(sigmas, deltas)):
        plt.plot(t, delta[1:], c[i], label=fr'$\sigma={s}$')

    plt.legend()
    plt.grid()
    plt.xlim(t[0], t[-1])
    #plt.xlabel(r'$t$')
    #plt.ylabel(r'$E$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()

def delta_different_tau(ts, deltas, taus, dir, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    
    for i, (tau, t, delta) in enumerate(zip(taus, ts, deltas)):
        plt.plot(t, delta[1:], c[i], label=fr'$\tau={tau}$')

    plt.legend()
    plt.grid()
    plt.xlim(t[0], t[-1])
    #plt.xlabel(r'$t$')
    #plt.ylabel(r'$E$')

    plt.savefig(join(dir, f'{filename}.png'))
    
    plt.close()


if __name__ == '__main__':
    pass
