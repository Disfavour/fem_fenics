from os.path import dirname, join
import matplotlib.pyplot as plt
from sw_1d import *
from sw_1d_smooth import *

l = ('-', '--', '-.', ':')
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
lw = [0.9, 1.2, 1.5]
dpi = 300

def t2_h_dif_t(x, h, dif_t, fname):
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

def t2_u_dif_t(x, u, dif_t, fname):
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

def t1_h_dif_t(x, h, he, dif_t, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, t, n, e in zip(colors, dif_t, h, he):
        plt.plot(x, n, f'-{c}', label=fr'$t={t}$')
        plt.plot(x, e, f':{c}')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$h$')
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def t1_u_dif_t(x, u, ue, dif_t, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, t, n, e in zip(colors, dif_t, u, ue):
        plt.plot(x, n, f'-{c}', label=fr'$t={t}$')
        plt.plot(x, e, f':{c}')
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
    taus = [0.01, 0.005, 0.0025]
    ms = [100, 200, 400]

    # for s in sigmas:
    #     for tau in taus:
    #         for m in ms:
    #             try:
    #                 x, h, u, he, ue, t, E = sw_1d(s=s, tau=tau, mesh_size=m, T=T1, dif_t=dif_t1, critical_time=ct)
    #             except Exception:
    #                 continue
    #             t1_h_dif_t(x, h, he, dif_t1, join(images_dir, f't1_h_dif_t_s{s}_tau{tau}_m{m}.png'))
    #             t1_h_dif_t(x, u, ue, dif_t1, join(images_dir, f't1_u_dif_t_s{s}_tau{tau}_m{m}.png'))

    #             x, h, u, t, E = sw_1d_smooth(s=s, tau=tau, mesh_size=m, T=T2, dif_t=dif_t2)
    #             t2_h_dif_t(x, h, dif_t2, join(images_dir, f't2_h_dif_t_s{s}_tau{tau}_m{m}.png'))
    #             t2_u_dif_t(x, u, dif_t2, join(images_dir, f't2_u_dif_t_s{s}_tau{tau}_m{m}.png'))

    # for s in sigmas:
    #     for tau in taus:
    #         ts, Es = [], []
    #         for m in ms:
    #             x, h, u, t, E = sw_1d_smooth(s=s, tau=tau, mesh_size=m, T=T2, dif_t=dif_t2)
    #             ts.append(t)
    #             Es.append(E)
    #         E_different_m(ts, Es, ms, join(images_dir, f't2_E_dif_m_s{s}_tau{tau}.png'))

    #         ts, Es = [], []
    #         msc = []
    #         for m in ms:
    #             try:
    #                 x, h, u, he, ue, t, E = sw_1d(s=s, tau=tau, mesh_size=m, T=T1, dif_t=dif_t1, critical_time=ct)
    #             except Exception:
    #                 continue
    #             ts.append(t)
    #             Es.append(E)
    #             msc.append(m)
    #         if msc:
    #             E_different_m(ts, Es, msc, join(images_dir, f't1_E_dif_m_s{s}_tau{tau}.png'))
    
    # for s in sigmas:
    #     for m in ms:
    #         ts, Es = [], []
    #         for tau in taus:
    #             x, h, u, t, E = sw_1d_smooth(s=s, tau=tau, mesh_size=m, T=T2, dif_t=dif_t2)
    #             ts.append(t)
    #             Es.append(E)
    #         E_different_tau(ts, Es, taus, join(images_dir, f't2_E_dif_tau_s{s}_m{m}.png'))

    #         ts, Es = [], []
    #         tausc = []
    #         for tau in taus:
    #             try:
    #                 x, h, u, he, ue, t, E = sw_1d(s=s, tau=tau, mesh_size=m, T=T1, dif_t=dif_t1, critical_time=ct)
    #             except Exception:
    #                 continue
    #             ts.append(t)
    #             Es.append(E)
    #             tausc.append(tau)
    #         if tausc:
    #             E_different_tau(ts, Es, tausc, join(images_dir, f't1_E_dif_tau_s{s}_m{m}.png'))
        
    for tau in taus:
        for m in ms:
            # ts, Es = [], []
            # for s in sigmas:
            #     x, h, u, t, E = sw_1d_smooth(s=s, tau=tau, mesh_size=m, T=T2, dif_t=dif_t2)
            #     ts.append(t)
            #     Es.append(E)
            # E_different_s(ts, Es, sigmas, join(images_dir, f't2_E_dif_s_tau{tau}_m{m}.png'))

            ts, Es = [], []
            sigmasc = []
            for s in sigmas:
                try:
                    x, h, u, he, ue, t, E = sw_1d(s=s, tau=tau, mesh_size=m, T=T1, dif_t=dif_t1, critical_time=ct)
                except Exception:
                    continue
                ts.append(t)
                Es.append(E)
                sigmasc.append(s)
            if sigmasc:
                E_different_s(ts, Es, sigmasc, join(images_dir, f't1_E_dif_s_tau{tau}_m{m}.png'))
