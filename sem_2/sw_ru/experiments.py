from multiprocessing import Pool
from os.path import join, dirname
from calcs.sw_r_u import sw_r_u
from calcs.sw_r_ru import sw_r_ru
from calcs.sw_s_u import sw_s_u
from calcs.sw_s_su import sw_s_su
from plotting import density_at_different_moments, E, E_compare
import numpy as np

images = join(dirname(__file__), 'images')
data = join(dirname(__file__), 'data')


# def r_u(ms, tau, degree_r, degree_u, sigma):
#     ts, Es, r_min, r_max, r_00, xs, r_ = sw_r_u(ms, tau, degree_r, degree_u, sigma)
#     density_at_different_moments(xs, *r_, join(images, f'density_r_u_ms{ms}_tau{tau}_d{degree_r}{degree_u}_sigma{sigma}.png'))
#     E(ts, Es, join(images, f'E_r_u_ms{ms}_tau{tau}_d{degree_r}{degree_u}_sigma{sigma}.png'))


# def r_ru(ms, tau, degree_r, degree_w, sigma):
#     ts, Es, r_min, r_max, r_00, xs, r_ = sw_r_ru(ms, tau, degree_r, degree_w, sigma)
#     density_at_different_moments(xs, *r_, join(images, f'density_r_ru_ms{ms}_tau{tau}_d{degree_r}{degree_w}_sigma{sigma}.png'))
#     E(ts, Es, join(images, f'E_r_ru_ms{ms}_tau{tau}_d{degree_r}{degree_w}_sigma{sigma}.png'))


# def s_u(ms, tau, degree_s, degree_u, sigma):
#     ts, Es, r_min, r_max, r_00, xs, r_ = sw_s_u(ms, tau, degree_s, degree_u, sigma)
#     density_at_different_moments(xs, *r_, join(images, f'density_s_u_ms{ms}_tau{tau}_d{degree_s}{degree_u}_sigma{sigma}.png'))
#     E(ts, Es, join(images, f'E_s_u_ms{ms}_tau{tau}_d{degree_s}{degree_u}_sigma{sigma}.png'))


# def s_su(ms, tau, degree_s, degree_w, sigma):
#     ts, Es, r_min, r_max, r_00, xs, r_ = sw_s_su(ms, tau, degree_s, degree_w, sigma)
#     density_at_different_moments(xs, *r_, join(images, f'density_s_su_ms{ms}_tau{tau}_d{degree_s}{degree_w}_sigma{sigma}.png'))
#     E(ts, Es, join(images, f'E_s_su_ms{ms}_tau{tau}_d{degree_s}{degree_w}_sigma{sigma}.png'))


def r_u(ms, tau, d1, d2, sigma):
    ts, Es, r_min, r_max, r_00, xs, r_ = sw_r_u(ms, tau, d1, d2, sigma)
    np.save(join(data, f'x_r_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((xs, *r_)))
    np.save(join(data, f't_r_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((ts, Es, r_min, r_max, r_00)))

def r_ru(ms, tau, d1, d2, sigma):
    ts, Es, r_min, r_max, r_00, xs, r_ = sw_r_ru(ms, tau, d1, d2, sigma)
    np.save(join(data, f'x_r_ru_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((xs, *r_)))
    np.save(join(data, f't_r_ru_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((ts, Es, r_min, r_max, r_00)))

def s_u(ms, tau, d1, d2, sigma):
    ts, Es, r_min, r_max, r_00, xs, r_ = sw_s_u(ms, tau, d1, d2, sigma)
    np.save(join(data, f'x_s_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((xs, *r_)))
    np.save(join(data, f't_s_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((ts, Es, r_min, r_max, r_00)))

def s_su(ms, tau, d1, d2, sigma):
    ts, Es, r_min, r_max, r_00, xs, r_ = sw_s_su(ms, tau, d1, d2, sigma)
    np.save(join(data, f'x_s_su_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((xs, *r_)))
    np.save(join(data, f't_s_su_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}'), np.array((ts, Es, r_min, r_max, r_00)))

def make_plots(ms, tau, d1, d2, sigma):
    x_r_u = np.load(join(data, f'x_r_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))
    t_r_u = np.load(join(data, f't_r_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))

    x_r_ru = np.load(join(data, f'x_r_ru_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))
    t_r_ru = np.load(join(data, f't_r_ru_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))

    x_s_u = np.load(join(data, f'x_s_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))
    t_s_u = np.load(join(data, f't_s_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))

    x_s_su = np.load(join(data, f'x_s_su_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))
    t_s_su = np.load(join(data, f't_s_su_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.npy'))

    E_compare(t_r_u[0], t_r_u[1], t_r_ru[1], t_s_u[1], t_s_su[1], 'r_u', 'r_ru', 's_u', 's_su', join(images, f'E_compare_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))

    density_at_different_moments(x_r_u[0], x_r_u[1], x_r_u[2], x_r_u[3], join(images, f'density_r_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))
    density_at_different_moments(x_r_ru[0], x_r_ru[1], x_r_ru[2], x_r_ru[3], join(images, f'density_r_ru_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))
    density_at_different_moments(x_s_u[0], x_s_u[1], x_s_u[2], x_s_u[3], join(images, f'density_s_u_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))
    density_at_different_moments(x_s_su[0], x_s_su[1], x_s_su[2], x_s_su[3], join(images, f'density_s_su_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))


def run():
    mesh_sizes = (100,)
    taus = (0.01,)
    degree1 = (2, 1)
    degree2 = (2, 1)
    sigmas = (0.5, 1)

    params = tuple((ms, tau, d1, d2, sigma) for ms in mesh_sizes for d2 in degree2 for d1 in degree1 for tau in taus for sigma in sigmas)

    #r_ru(50, 0.02, 1, 1, 1)

    # for p in params:
    #     r_u(*p)

    with Pool() as p:
        #print(p._processes)
        res1 = p.starmap_async(r_u, params)
        res2 = p.starmap_async(r_ru, params)
        res3 = p.starmap_async(s_u, params)
        res4 = p.starmap_async(s_su, params)
        res1.wait()
        res2.wait()
        res3.wait()
        res4.wait()

        res5 = p.starmap_async(make_plots, params)
        res5.wait()
    
    # for i in params:
    #     make_plots(*i)


def special():
    ms = 100
    tau = 0.01
    d1 = d2 = 1
    sigma = 0.5
    ts, Es, r_min, r_max, r_00, xs, r_ = sw_s_su(100, 0.01, 1, 1, 0.5)
    density_at_different_moments(xs, *r_, join(images, f'density_s_su_fix_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))
    E(ts, Es, join(images, f'E_s_su_fix_ms{ms}_tau{tau}_d{d1}{d2}_sigma{sigma}.png'))


if __name__ == '__main__':
    #run()
    special()