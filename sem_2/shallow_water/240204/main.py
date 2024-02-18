import fem_with_analytic
import triple_layer
import triple_layer_improved
import plotting
from os.path import dirname, join, exists
import numpy as np
from multiprocessing import Pool
import time


def store_data(hl, hr, T, mesh_size, tau, theta, sigma, ts2store, get_solution, data_dir):
    try:
        results = get_solution(hl, hr, T, mesh_size, tau, theta, sigma, ts2store)
    except:
        return
    t, m, E, m_e, E_e, err_h, err_u, x, h, u, x_e, h_e, u_e = results
    np.savez_compressed(join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz'),
                        t=t, m=m, E=E, m_e=m_e, E_e=E_e, err_h=err_h, err_u=err_u,
                        x=x, h=h, u=u, x_e=x_e, h_e=h_e, u_e=u_e,
                        ts=np.array(ts2store))


def plot_dif_times(hl, hr, mesh_size, tau, theta, sigma, data_dir, images_dir):
    fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz')
    if exists(fname):
        npz = np.load(fname)
        plotting.dif_times(npz['x'], npz['h'], npz['x_e'], npz['h_e'], r'$h$', [f'$t={i}$' for i in npz['ts']],
                           join(images_dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.png'))
        plotting.dif_times(npz['x'], npz['u'], npz['x_e'], npz['u_e'], r'$u$', [f'$t={i}$' for i in npz['ts']],
                           join(images_dir, f'hl{hl}_hr{hr}_u_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.png'))
    

def plot_dif_ms(hl, hr, mesh_sizes, tau, theta, sigma, data_dir, images_dir):
    d, ms_ = [], []
    for mesh_size in reversed(mesh_sizes):
        fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz')
        if exists(fname):
            d.append(np.load(fname))
            ms_.append(mesh_size)
    if len(d) > 0:
        plotting.base([i['t'] for i in d], [i['err_h'] for i in d], r'$t$', r'$||h_{exact} - h_{numerical}||_2$', [rf'$M={i}$' for i in ms_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d], [i['err_u'] for i in d], r'$t$', r'$||u_{exact} - u_{numerical}||_2$', [rf'$M={i}$' for i in ms_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d] + [d[0]['t']], [i['E'] for i in d] + [d[0]['E_e']], r'$t$', r'$E$', [rf'$M={i}$' for i in ms_] + ['exact'],
                      join(images_dir, f'hl{hl}_hr{hr}_E_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png'))


def plot_dif_tau(hl, hr, mesh_size, taus, theta, sigma, data_dir, images_dir):
    d, tau_ = [], []
    for tau in taus:
        fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz')
        if exists(fname):
            d.append(np.load(fname))
            tau_.append(tau)
    if len(d) > 0:
        plotting.base([i['t'] for i in d], [i['err_h'] for i in d], r'$t$', r'$||h_{exact} - h_{numerical}||_2$', [rf'$\tau={i}$' for i in tau_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_h_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d], [i['err_u'] for i in d], r'$t$', r'$||u_{exact} - u_{numerical}||_2$', [rf'$\tau={i}$' for i in tau_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_u_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d] + [d[0]['t']], [i['E'] for i in d] + [d[0]['E_e']], r'$t$', r'$E$', [rf'$\tau={i}$' for i in tau_] + ['exact'],
                      join(images_dir, f'hl{hl}_hr{hr}_E_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png'))


def plot_dif_theta(hl, hr, mesh_size, tau, thetas, sigma, data_dir, images_dir):
    d, theta_ = [], []
    for theta in thetas:
        fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz')
        if exists(fname):
            d.append(np.load(fname))
            theta_.append(theta)
    if len(d) > 0:
        plotting.base([i['t'] for i in d], [i['err_h'] for i in d], r'$t$', r'$||h_{exact} - h_{numerical}||_2$', [rf'$\theta={i}$' for i in theta_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_h_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d], [i['err_u'] for i in d], r'$t$', r'$||u_{exact} - u_{numerical}||_2$', [rf'$\theta={i}$' for i in theta_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_u_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png'))
        plotting.base([i['t'] for i in d] + [d[0]['t']], [i['E'] for i in d] + [d[0]['E_e']], r'$t$', r'$E$', [rf'$\theta={i}$' for i in theta_] + ['exact'],
                      join(images_dir, f'hl{hl}_hr{hr}_E_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png'))


def plot_dif_sigma(hl, hr, mesh_size, tau, theta, sigmas, data_dir, images_dir):
    d, s_ = [], []
    for sigma in sigmas:
        fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.npz')
        if exists(fname):
            d.append(np.load(fname))
            s_.append(sigma)
    if len(d) > 0:
        plotting.base([i['t'] for i in d], [i['err_h'] for i in d], r'$t$', r'$||h_{exact} - h_{numerical}||_2$', [rf'$\sigma={i}$' for i in s_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_h_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png'))
        plotting.base([i['t'] for i in d], [i['err_u'] for i in d], r'$t$', r'$||u_{exact} - u_{numerical}||_2$', [rf'$\sigma={i}$' for i in s_],
                      join(images_dir, f'hl{hl}_hr{hr}_err_u_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png'))
        plotting.base([i['t'] for i in d] + [d[0]['t']], [i['E'] for i in d] + [d[0]['E_e']], r'$t$', r'$E$', [rf'$\sigma={i}$' for i in s_] + ['exact'],
                      join(images_dir, f'hl{hl}_hr{hr}_E_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png'))


def experiments(hl, hr, T, mesh_sizes, taus, thetas, sigmas, get_solution, data_dir, images_dir):
    ts = np.arange(0, T+taus[-1]/2, taus[-1])
    ts2store = [ts[2], ts[ts.shape[0]//2+1], ts[-1]]

    params = [[hl, hr, T, mesh_size, tau, theta, sigma, ts2store, get_solution, data_dir]
              for mesh_size in mesh_sizes for tau in taus for theta in thetas for sigma in sigmas]
    
    with Pool() as p:
        res1 = p.starmap_async(store_data, params)
        res1.wait()

    params_dif = [[hl, hr, mesh_size, tau, theta, sigma, data_dir, images_dir]
              for mesh_size in mesh_sizes for tau in taus for theta in thetas for sigma in sigmas]

    params_dif_ms = [[hl, hr, mesh_sizes, tau, theta, sigma, data_dir, images_dir]
              for tau in taus for theta in thetas for sigma in sigmas]
    params_dif_tau = [[hl, hr, mesh_size, taus, theta, sigma, data_dir, images_dir]
              for mesh_size in mesh_sizes for theta in thetas for sigma in sigmas]
    params_dif_theta = [[hl, hr, mesh_size, tau, thetas, sigma, data_dir, images_dir]
              for mesh_size in mesh_sizes for tau in taus for sigma in sigmas]
    params_dif_sigma = [[hl, hr, mesh_size, tau, theta, sigmas, data_dir, images_dir]
              for mesh_size in mesh_sizes for tau in taus for theta in thetas]
    
    with Pool() as p:
        tmp = []
        tmp.append(p.starmap_async(plot_dif_times, params_dif))
        tmp.append(p.starmap_async(plot_dif_ms, params_dif_ms))
        tmp.append(p.starmap_async(plot_dif_tau, params_dif_tau))
        tmp.append(p.starmap_async(plot_dif_theta, params_dif_theta))
        tmp.append(p.starmap_async(plot_dif_sigma, params_dif_sigma))
        list(map(lambda x: x.wait(), tmp))


def test_1(hl, hr, T, base_data_dir, base_images_dir):
    mesh_sizes = [3200, 1600, 800, 400, 200, 100]
    taus = [0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08]
    sigmas = [1.0]
    thetas = [0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]

    triple_layer_data = join(base_data_dir, 'triple_layer')
    triple_layer_images = join(base_images_dir, 'triple_layer')

    experiments(hl, hr, T, mesh_sizes, taus, thetas, sigmas, triple_layer.get_solution, triple_layer_data, triple_layer_images)


if __name__ == '__main__':
    base_dir = dirname(__file__)
    data_dir = join(base_dir, 'data')
    images_dir = join(base_dir, 'images')
    moments = []

    name = 'triple_layer'
    mesh_sizes = [3200, 1600, 800, 400, 200, 100]
    taus = [0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08]
    thetas = [0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    sigmas = [1.0]

    moments.append(time.time())
    #experiments(10, 1, 0.4, mesh_sizes, taus, thetas, sigmas, triple_layer.get_solution, join(data_dir, name), join(images_dir, name))
    moments.append(time.time())
    #experiments(2, 1, 0.9, mesh_sizes, taus, thetas, sigmas, triple_layer.get_solution, join(data_dir, name), join(images_dir, name))
    moments.append(time.time())

    name = 'triple_layer_improved'
    mesh_sizes = [3200, 1600, 800, 400, 200, 100]
    taus = [0.00125, 0.0025, 0.005, 0.01, 0.02, 0.04, 0.08]
    thetas = [0.5]
    sigmas = [0.125, 0.167, 0.25, 0.5, 1.0]

    experiments(10, 1, 0.4, mesh_sizes, taus, thetas, sigmas, triple_layer_improved.get_solution, join(data_dir, name), join(images_dir, name))
    moments.append(time.time())
    #experiments(2, 1, 0.9, mesh_sizes, taus, thetas, sigmas, triple_layer_improved.get_solution, join(data_dir, name), join(images_dir, name))
    moments.append(time.time())

    moments = np.array(moments)
    print(moments)
    print(moments[1:] - moments[:-1])

    # вернуть сетку 100к?


    #t, m, E, m_e, E_e, err_h, err_u, x, h, u, x_e, h_e, u_e = fem_with_analytic.get_solution(10, 1, 0.4, 200, 0.005, 1.0, [0.1, 0.5, 0.9])
    #plotting.dif_times(x, h, x_e, h_e, r'$h$', [f'$t={i}$' for i in [0.1, 0.5, 0.9]], join(images_dir, f'test.png'))
    #plotting.base([t, t], [m, m_e], r'$t$', r'$m$', ['m', 'm_e'], join(images_dir, f'test.png'))
    