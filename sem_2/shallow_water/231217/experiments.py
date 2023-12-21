from os.path import dirname, join
import matplotlib.pyplot as plt
from sw_2d_v1 import *
from sw_2d_v2 import *
from sw_2d_v3 import *

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

def hp_dif_t(x, hs, hps, dif_t, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for c, t, h, hp in zip(colors, dif_t, hs, hps):
        plt.plot(x, abs(h - hp), c, label=fr'$t={t}$')
        #plt.plot(x, hp, f'--{c}')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$|h - h(p)|$')
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_m(ts, ms, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(ts, ms)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$m$')
    plt.xlim(ts[0], ts[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()

def plot_E(ts, Es, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    plt.plot(ts, Es)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.xlim(ts[0], ts[-1])
    plt.grid()
    plt.savefig(fname)
    plt.close()

if __name__ == '__main__':
    base_dir = dirname(__file__)
    images_dir = join(base_dir, 'images')

    sigmas = [0.5, 0.75, 1, 1.5]
    taus = [0.01, 0.005, 0.0025, 0.00125]
    ms = [100, 200, 400, 800]

    T = 2
    dif_t = (0.1, 0.7, 1.3, 1.9)

    sigmas = [0.5, 0.75, 1, 1.5]
    taus = [0.005]
    mesh_sizes = [200]

    # for s in sigmas:
    #     for tau in taus:
    #         for mesh_size in mesh_sizes:
    #             ts, ms, Es, x, hs, hsp = sw_2d_v3(s, tau, mesh_size, T, dif_t)
    #             for t, h, hp in zip(dif_t, hs, hsp):
    #                 hp_dif_t(x, [h], [hp], [t], join(images_dir, f'v{3}_hp_{t}_s{s}_tau{tau}_m{mesh_size}.png'))

    for s in sigmas:
        for tau in taus:
            for mesh_size in mesh_sizes:
                for v, sw in zip([1, 2, 3], [sw_2d_v1, sw_2d_v2, sw_2d_v3]):
                #for v, sw in zip([3], [sw_2d_v3]):
                    try:
                        ts, ms, Es, x, hs, hsp = sw(s, tau, mesh_size, T, dif_t)
                    except Exception:
                        continue
                    h_dif_t(x, hs, dif_t, join(images_dir, f'v{v}_h_s{s}_tau{tau}_m{mesh_size}.png'))
                    plot_m(ts, ms, join(images_dir, f'v{v}_m_s{s}_tau{tau}_m{mesh_size}.png'))
                    plot_E(ts, Es, join(images_dir, f'v{v}_E_s{s}_tau{tau}_m{mesh_size}.png'))
