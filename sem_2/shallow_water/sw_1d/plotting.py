import matplotlib.pyplot as plt
import numpy as np
from os.path import join


l = ('-', '--', ':') # '-.'
c = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
lw = [1.5, 1.0, 0.5]
dpi = 300

def h_at_different_moments(x, h1, h2, h3, time_moments, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(x[2], h3[0], ':r', lw=lw[0])
    plt.plot(x[2], h3[1], ':g', lw=lw[0])
    plt.plot(x[2], h3[2], ':b', lw=lw[0])

    plt.plot(x[1], h2[0], '--r', lw=lw[1])
    plt.plot(x[1], h2[1], '--g', lw=lw[1])
    plt.plot(x[1], h2[2], '--b', lw=lw[1])

    plt.plot(x[0], h1[0], '-r', lw=lw[2], label=fr'$t={time_moments[0]}$')
    plt.plot(x[0], h1[1], '-g', lw=lw[2], label=fr'$t={time_moments[1]}$')
    plt.plot(x[0], h1[2], '-b', lw=lw[2], label=fr'$t={time_moments[2]}$')

    plt.xlabel(r'$x$')
    plt.ylabel(r'$h$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()

def u_at_different_moments(x, u1, u2, u3, time_moments, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(x[2], u3[0], ':r', lw=lw[0])
    plt.plot(x[2], u3[1], ':g', lw=lw[0])
    plt.plot(x[2], u3[2], ':b', lw=lw[0])

    plt.plot(x[1], u2[0], '--r', lw=lw[1])
    plt.plot(x[1], u2[1], '--g', lw=lw[1])
    plt.plot(x[1], u2[2], '--b', lw=lw[1])

    plt.plot(x[0], u1[0], '-r', lw=lw[2], label=fr'$t={time_moments[0]}$')
    plt.plot(x[0], u1[1], '-g', lw=lw[2], label=fr'$t={time_moments[1]}$')
    plt.plot(x[0], u1[2], '-b', lw=lw[2], label=fr'$t={time_moments[2]}$')
    
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()


def E_different_tau(ts, Es, taus, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], Es[0], lw=lw[0], label=fr'$t={taus[0]}$')
    plt.plot(ts[1], Es[1], lw=lw[1], label=fr'$t={taus[1]}$')
    plt.plot(ts[2], Es[2], lw=lw[2], label=fr'$t={taus[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()

def E_different_sigma(ts, Es, sigmas, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], Es[0], lw=lw[0], label=fr'$\sigma={sigmas[0]}$')
    plt.plot(ts[1], Es[1], lw=lw[1], label=fr'$\sigma={sigmas[1]}$')
    plt.plot(ts[2], Es[2], lw=lw[2], label=fr'$\sigma={sigmas[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()

def E_different_ms(ts, Es, mesh_sizes, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], Es[0], lw=lw[0], label=fr'$mesh \ size={mesh_sizes[0]}$')
    plt.plot(ts[1], Es[1], lw=lw[1], label=fr'$mesh \ size={mesh_sizes[1]}$')
    plt.plot(ts[2], Es[2], lw=lw[2], label=fr'$mesh \ size={mesh_sizes[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()


def courant_different_tau(ts, courant, taus, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], courant[0], lw=lw[0], label=fr'$t={taus[0]}$')
    plt.plot(ts[1], courant[1], lw=lw[1], label=fr'$t={taus[1]}$')
    plt.plot(ts[2], courant[2], lw=lw[2], label=fr'$t={taus[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$courant$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()

def courant_different_sigma(ts, courant, sigmas, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], courant[0], lw=lw[0], label=fr'$\sigma={sigmas[0]}$')
    plt.plot(ts[1], courant[1], lw=lw[1], label=fr'$\sigma={sigmas[1]}$')
    plt.plot(ts[2], courant[2], lw=lw[2], label=fr'$\sigma={sigmas[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$courant$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()

def courant_different_ms(ts, courant, mesh_sizes, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(ts[0], courant[0], lw=lw[0], label=fr'$mesh \ size={mesh_sizes[0]}$')
    plt.plot(ts[1], courant[1], lw=lw[1], label=fr'$mesh \ size={mesh_sizes[1]}$')
    plt.plot(ts[2], courant[2], lw=lw[2], label=fr'$mesh \ size={mesh_sizes[2]}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$courant$')
    plt.grid()
    plt.legend()
    plt.savefig(fname)
    plt.close()