from os.path import dirname, join
import matplotlib.pyplot as plt
from sw_1d import *
from sw_1d_smooth import *

l = ('-', '--', '-.', ':')
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
lw = [0.9, 1.2, 1.5]
dpi = 300

def h_dif_t(x, h, dif_t, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, t, n in zip(colors, dif_t, h):
        plt.plot(x, n, c, label=fr'$t={t}$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$h$')
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def u_dif_t(x, u, dif_t, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, t, n in zip(colors, dif_t, u):
        plt.plot(x, n, c, label=fr'$t={t}$')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def E_different_s(ts, Es, sigmas, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, s, t, E in zip(colors, sigmas, ts, Es):
        plt.plot(t, E, c, label=fr'$\sigma={s}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def E_different_tau(ts, Es, taus, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, tau, t, E in zip(colors, taus, ts, Es):
        plt.plot(t, E, c, label=fr'$\tau={tau}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def E_different_m(ts, Es, ms, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, m, t, E in zip(colors, ms, ts, Es):
        plt.plot(t, E, c, label=fr'$m={m}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()


if __name__ == '__main__':
    base_dir = dirname(__file__)
    images_dir = join(base_dir, 'images')

    # sigma
    sigmas = [0.75, 1, 1.5]#[0.75, 1, 1.5, 2]
    tau = 0.005
    mesh_size = 200
    T1 = 1
    T2 = 2
    
    ct = 0.5 # подобрано

    dif_t1 = [0.1, 0.4, 0.7, 1]
    dif_t2 = [0.1, 0.7, 1.3, 1.9]

    sigmas = [0.5, 0.75, 1, 1.5]
    taus = [0.01, 0.005, 0.0025, 0.00125]
    ms = [100, 200, 400, 800]

    for s in sigmas:
        for tau in taus:
            for m in ms:
                for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
                    for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
                        try:
                            x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
                        except Exception:
                            continue
                        
                        h_dif_t(x, h, dif_t, join(images_dir, f't{tn}_hl{hl}_h_dif_t_s{s}_tau{tau}_m{m}.png'))
                        u_dif_t(x, u, dif_t, join(images_dir, f't{tn}_hl{hl}_u_dif_t_s{s}_tau{tau}_m{m}.png'))

    for s in sigmas:
        for tau in taus:
            for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
                for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
                    ts, Es, msc = [], [], []
                    for m in ms:
                        try:
                            x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
                        except Exception:
                            continue
                        ts.append(t)
                        Es.append(E)
                        msc.append(m)
                    if msc:
                        E_different_m(ts, Es, msc, join(images_dir, f't{tn}_hl{hl}_E_dif_m_s{s}_tau{tau}.png'))
    
    for s in sigmas:
        for m in ms:
            for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
                for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
                    ts, Es, tausc = [], [], []
                    for tau in taus:
                        try:
                            x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
                        except Exception:
                            continue
                        ts.append(t)
                        Es.append(E)
                        tausc.append(tau)
                    if tausc:
                        E_different_tau(ts, Es, tausc, join(images_dir, f't{tn}_hl{hl}_E_dif_tau_s{s}_m{m}.png'))

    for tau in taus:
        for m in ms:
            for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
                for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
                    ts, Es, sigmasc = [], [], []
                    for s in sigmas:
                        try:
                            x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
                        except Exception:
                            continue
                        ts.append(t)
                        Es.append(E)
                        sigmasc.append(s)
                    if sigmasc:
                        E_different_s(ts, Es, sigmasc, join(images_dir, f't{tn}_hl{hl}_E_dif_s_tau{tau}_m{m}.png'))
