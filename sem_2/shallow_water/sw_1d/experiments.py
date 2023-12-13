from os.path import dirname, join
from sw_1d import sw_1d
from plotting import *

base_dir = dirname(__file__)
images_dir = join(base_dir, 'images')


def run_test_problem(sigmas, taus, mesh_sizes, T, time_moments, depths):
    for depth in depths:
        for mesh_size in mesh_sizes:
            for sigma in sigmas:
                xs, hs, us, ts, Es, courants = [], [], [], [], [], []
                for tau in taus:
                    x, h, u, time_moments, t, E, courant = sw_1d(s=sigma, tau=tau, mesh_size=mesh_size, T=T, time_moments=time_moments, depth=depth)
                    xs.append(x)
                    hs.append(h)
                    us.append(u)
                    ts.append(t)
                    Es.append(E)
                    courants.append(courant)
                h_at_different_moments(xs, *hs, time_moments, join(images_dir, f'depth{depth}_h_at_different_moments_different_tau_s{sigma}_ms{mesh_size}.png'))
                u_at_different_moments(xs, *us, time_moments, join(images_dir, f'depth{depth}_u_at_different_moments_different_tau_s{sigma}_ms{mesh_size}.png'))
                E_different_tau(ts, Es, taus, join(images_dir, f'depth{depth}_E_different_tau_s{sigma}_ms{mesh_size}.png'))
                courant_different_tau(ts, courants, taus, join(images_dir, f'depth{depth}_courant_different_tau_s{sigma}_ms{mesh_size}.png'))
        
        for mesh_size in mesh_sizes:
            for tau in taus:
                xs, hs, us, ts, Es, courants = [], [], [], [], [], []
                for sigma in sigmas:
                    x, h, u, time_moments, t, E, courant = sw_1d(s=sigma, tau=tau, mesh_size=mesh_size, T=T, time_moments=time_moments, depth=depth)
                    xs.append(x)
                    hs.append(h)
                    us.append(u)
                    ts.append(t)
                    Es.append(E)
                    courants.append(courant)
                h_at_different_moments(xs, *hs, time_moments, join(images_dir, f'depth{depth}_h_at_different_moments_different_sigma_tau{tau}_ms{mesh_size}.png'))
                u_at_different_moments(xs, *us, time_moments, join(images_dir, f'depth{depth}_u_at_different_moments_different_sigma_tau{tau}_ms{mesh_size}.png'))
                E_different_sigma(ts, Es, sigmas, join(images_dir, f'depth{depth}_E_different_sigma_tau{tau}_ms{mesh_size}.png'))
                courant_different_sigma(ts, courants, sigmas, join(images_dir, f'depth{depth}_courant_different_sigma_tau{tau}_ms{mesh_size}.png'))

        for sigma in sigmas:
            for tau in taus:
                xs, hs, us, ts, Es, courants = [], [], [], [], [], []
                for mesh_size in mesh_sizes:
                    x, h, u, time_moments, t, E, courant = sw_1d(s=sigma, tau=tau, mesh_size=mesh_size, T=T, time_moments=time_moments, depth=depth)
                    xs.append(x)
                    hs.append(h)
                    us.append(u)
                    ts.append(t)
                    Es.append(E)
                    courants.append(courant)
                h_at_different_moments(xs, *hs, time_moments, join(images_dir, f'depth{depth}_h_at_different_moments_different_ms_s{sigma}_tau{tau}.png'))
                u_at_different_moments(xs, *us, time_moments, join(images_dir, f'depth{depth}_u_at_different_moments_different_ms_s{sigma}_tau{tau}.png'))
                E_different_ms(ts, Es, mesh_sizes, join(images_dir, f'depth{depth}_E_different_ms_s{sigma}_tau{tau}.png'))
                courant_different_ms(ts, courants, mesh_sizes, join(images_dir, f'depth{depth}_courant_different_ms_s{sigma}_tau{tau}.png'))

if __name__ == '__main__':
    sigma = (0.75, 1, 1.5)
    taus = (0.0025, 0.005, 0.01)
    mesh_size = (400, 200, 100)
    #mesh_size = (800, 400, 200)
    T = 2
    time_moments = (0.5, 1.0, 1.5)
    depths = (0.0625,)
    #run_test_problem(sigma, tau, mesh_size, T, time_moments, depths)

    for tau, ms in zip(taus, mesh_size):
        x, h, u, time_moments, t, E, courant, delta = sw_1d(s=1, tau=tau, mesh_size=ms, T=5, time_moments=[], depth=1)
        plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
        plt.plot(t, E)
        plt.xlabel(r'$t$')
        plt.ylabel(r'$E$')
        plt.grid()
        plt.savefig(join(images_dir, f'E_s{1}_tau{tau}_ms{ms}.png'))
        plt.close()
        plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
        plt.plot(t, delta[1:])
        plt.xlabel(r'$t$')
        plt.ylabel(r'$\sum_{k=0}^n \tau \delta_h^k$')
        plt.grid()
        plt.savefig(join(images_dir, f'delta_s{1}_tau{tau}_ms{ms}.png'))
        plt.close()
