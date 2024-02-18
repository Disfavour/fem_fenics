from os.path import dirname, join
import matplotlib.pyplot as plt
import numpy as np
import v1, v2, v3, v4

l = ('-', '--', '-.', ':')
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
lw = [0.9, 1.2, 1.5]
dpi = 300

def plot_comp(t, comp, comp_ex, fname, ylbl):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, comp)
    plt.plot(t, comp_ex)
    plt.xlabel(r'$t$')
    plt.ylabel(ylbl)
    plt.xlim(t[0], t[-1])
    plt.legend(('numerical', 'exact'))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_dif(t, comp, comp_ex, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, np.abs(comp - comp_ex))
    plt.xlabel(r'$t$')
    plt.ylabel(f'$|numerical - exact|$')
    plt.xlim(t[0], t[-1])
    #plt.legend(('numerical', 'exact'))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_components(t, components, components_exact, fname):
    labels = [
        r'$\frac {h^{n+1} - h^n} {\tau} h^{n+\sigma}$',
        r'$\frac {\partial} {\partial x} (h u)  h^{n+\sigma}$',
        r'$\frac {h^{n+1}u^{n+1} - h^n u^n} {\tau}  u^{n+\sigma}$',
        r'$\frac {\partial} {\partial x} (h u u) u^{n+\sigma}$',
        r'$g h \frac {\partial h} {\partial x} u^{n+\sigma}$',
    ]
    
    for i, (label, comp, comp_ex) in enumerate(zip(labels, components, components_exact)):
        if i == 0:
            plot_comp(t[1:], comp[1:], comp_ex[1:], f'{fname}{i}.png', label)
            plot_dif(t[1:], comp[1:], comp_ex[1:], f'{fname}{i}_dif.png')
            continue
        plot_comp(t, comp, comp_ex, f'{fname}{i}.png', label)
        plot_dif(t, comp, comp_ex, f'{fname}{i}_dif.png')

def plot_m(t, m, me, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(t, m)
    plt.plot(t, me)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$m$')
    plt.xlim(t[0], t[-1])
    plt.legend(('numerical', 'exact'))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_E(t, E, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for E_ in E:
        #print(t.shape, E_.shape, len(E))
        plt.plot(t, E_)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(t[0], t[-1])
    plt.legend(range(1, len(E)+1))
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_dif_E(t, E, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for E_ in E[1:]:
        plt.plot(t, np.abs(E_ - E[0]))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$|E - E_1|$')
    plt.xlim(t[0], t[-1])
    plt.legend(range(2, len(E)+1))
    plt.grid()
    plt.savefig(fname)
    plt.close()

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

    sigmas = [0.5, 0.75, 1, 1.5]
    taus = [0.01, 0.005, 0.0025]
    mesh_sizes = [100, 200, 400]

    for ms in mesh_sizes:
        for tau in taus:
            for s in sigmas:
                Es = []
                try:
                    for sw_1d in (v1.sw_1d, v2.sw_1d, v3.sw_1d, v4.sw_1d):
                        t, m, E = sw_1d(ms, tau, s)
                        Es.append(E)
                except Exception:
                    continue
                plot_E(t, Es, join(images_dir, f'240107_ms{ms}_tau{tau}_s{s}_E.png'))
                plot_dif_E(t, Es, join(images_dir, f'240107_ms{ms}_tau{tau}_s{s}_dif_E.png'))

    # for s in sigmas:
    #     for tau in taus:
    #         for m in ms:
    #             for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
    #                 for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
    #                     try:
    #                         x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
    #                     except Exception:
    #                         continue
                        
    #                     h_dif_t(x, h, dif_t, join(images_dir, f't{tn}_hl{hl}_h_dif_t_s{s}_tau{tau}_m{m}.png'))
    #                     u_dif_t(x, u, dif_t, join(images_dir, f't{tn}_hl{hl}_u_dif_t_s{s}_tau{tau}_m{m}.png'))

    # for s in sigmas:
    #     for tau in taus:
    #         for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
    #             for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
    #                 ts, Es, msc = [], [], []
    #                 for m in ms:
    #                     try:
    #                         x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
    #                     except Exception:
    #                         continue
    #                     ts.append(t)
    #                     Es.append(E)
    #                     msc.append(m)
    #                 if msc:
    #                     E_different_m(ts, Es, msc, join(images_dir, f't{tn}_hl{hl}_E_dif_m_s{s}_tau{tau}.png'))
    
    # for s in sigmas:
    #     for m in ms:
    #         for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
    #             for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
    #                 ts, Es, tausc = [], [], []
    #                 for tau in taus:
    #                     try:
    #                         x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
    #                     except Exception:
    #                         continue
    #                     ts.append(t)
    #                     Es.append(E)
    #                     tausc.append(tau)
    #                 if tausc:
    #                     E_different_tau(ts, Es, tausc, join(images_dir, f't{tn}_hl{hl}_E_dif_tau_s{s}_m{m}.png'))

    # for tau in taus:
    #     for m in ms:
    #         for hl, dif_t, T in zip((3, 10), (dif_t2, dif_t1), (T2, T1)):
    #             for tn, sw in zip((1, 2), (sw_1d, sw_1d_smooth)):
    #                 ts, Es, sigmasc = [], [], []
    #                 for s in sigmas:
    #                     try:
    #                         x, h, u, t, E = sw(hl=hl, s=s, tau=tau, mesh_size=m, T=T, dif_t=dif_t)
    #                     except Exception:
    #                         continue
    #                     ts.append(t)
    #                     Es.append(E)
    #                     sigmasc.append(s)
    #                 if sigmasc:
    #                     E_different_s(ts, Es, sigmasc, join(images_dir, f't{tn}_hl{hl}_E_dif_s_tau{tau}_m{m}.png'))
