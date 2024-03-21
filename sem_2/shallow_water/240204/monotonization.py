import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np
from multiprocessing import Pool
import monotonization_1
import monotonization_2
import monotonization_all


def plot_dif_times(hl, hr, mesh_size, tau, alfa, gamma, data_dir, images_dir):
    fname = join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.npz')
    if exists(fname):
        npz = np.load(fname)
        plotting.dif_times(npz['x'], npz['h'], npz['x_e'], npz['h_e'], r'$h$', [f'$t={i}$' for i in npz['ts']],
                           join(images_dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.png'))
        plotting.dif_times(npz['x'], npz['u'], npz['x_e'], npz['u_e'], r'$u$', [f'$t={i}$' for i in npz['ts']],
                           join(images_dir, f'hl{hl}_hr{hr}_u_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.png'))
    
        # plotting.base([npz['t']], [npz['err_h']], r'$t$', r'$\varepsilon_h$', [],
        #               join(images_dir, f'hl{hl}_hr{hr}_err_h_tau{tau}_alfa{alfa}_gamma{gamma}.png'))
        # plotting.base([npz['t']], [npz['err_u']], r'$t$', r'$\varepsilon_u$', [],
        #               join(images_dir, f'hl{hl}_hr{hr}_err_u_tau{tau}_alfa{alfa}_gamma{gamma}.png'))


def store_data(hl, hr, T, mesh_size, tau, alfa, gamma, ts2store, calculate, data_dir):
    try:
        results = calculate(hl, hr, T, mesh_size, tau, alfa, gamma, ts2store)
    except:
        return
    t, m, E, m_e, E_e, err_h, err_u, x, h, u, x_e, h_e, u_e = results
    np.savez_compressed(join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.npz'),
                        t=t, m=m, E=E, m_e=m_e, E_e=E_e, err_h=err_h, err_u=err_u,
                        x=x, h=h, u=u, x_e=x_e, h_e=h_e, u_e=u_e,
                        ts=np.array(ts2store))


def experiments(hl, hr, T, mesh_sizes, taus, alfas, gammas, calculate, data_dir, images_dir):
    ts = np.arange(0, T+taus[-1]/2, taus[-1])

    ts2store = [ts[ts.shape[0]//3], ts[2*ts.shape[0]//3], ts[-1]]

    params = [[hl, hr, T, mesh_size, tau, alfa, gamma, ts2store, calculate, data_dir]
              for mesh_size in mesh_sizes for tau in taus for alfa in alfas for gamma in gammas
              if not exists(join(data_dir, f'hl{hl}_hr{hr}_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.npz'))]
    
    with Pool() as p:
        res1 = p.starmap_async(store_data, params)
        res1.wait()

    params_dif = ((hl, hr, mesh_size, tau, alfa, gamma, data_dir, images_dir)
              for mesh_size in mesh_sizes for tau in taus for alfa in alfas for gamma in gammas)

    # params_dif_ms = ((hl, hr, mesh_sizes, tau, theta, sigma, data_dir, images_dir)
    #           for tau in taus for theta in thetas for sigma in sigmas)
    # params_dif_tau = ((hl, hr, mesh_size, taus, theta, sigma, data_dir, images_dir)
    #           for mesh_size in mesh_sizes for theta in thetas for sigma in sigmas)
    # params_dif_theta = ((hl, hr, mesh_size, tau, thetas, sigma, data_dir, images_dir)
    #           for mesh_size in mesh_sizes for tau in taus for sigma in sigmas) if len(thetas) > 1 else ()
    # params_dif_sigma = ((hl, hr, mesh_size, tau, theta, sigmas, data_dir, images_dir)
    #           for mesh_size in mesh_sizes for tau in taus for theta in thetas) if len(sigmas) > 1 else ()
    
    with Pool() as p:
        tmp = []
        tmp.append(p.starmap_async(plot_dif_times, params_dif))
        # tmp.append(p.starmap_async(plot_dif_ms, params_dif_ms))
        # tmp.append(p.starmap_async(plot_dif_tau, params_dif_tau))
        # tmp.append(p.starmap_async(plot_dif_theta, params_dif_theta))
        # tmp.append(p.starmap_async(plot_dif_sigma, params_dif_sigma))
        list(map(lambda x: x.wait(), tmp))


def monotonization1(data_dir, images_dir):
    name = 'monotonization_1'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    alfas = [-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)]
    gammas = list(np.arange(-5, 5.1, 0.5))
    
    experiments(2, 1, 0.9, mesh_sizes, taus, alfas, gammas, monotonization_1.calculate, data_dir, images_dir)

def monotonization2(data_dir, images_dir):
    name = 'monotonization_2'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    alfas = [-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)]
    gammas = list(np.arange(-5, 5.1, 0.5))
    
    experiments(2, 1, 0.9, mesh_sizes, taus, alfas, gammas, monotonization_2.calculate, data_dir, images_dir)


def monotonization_full(data_dir, images_dir):
    name = 'monotonization_full'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    degrees = list(range(1, 4))
    T = 7.0
    sigma = 0.5
    alfas = [-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)]
    gammas = list(np.arange(-1, 4.1, 0.5))
    kappas = list(np.arange(0, 2.01, 0.25))
    ntests = (1, 2)


if __name__ == '__main__':
    import time
    base_dir = dirname(__file__)
    data_dir = join(base_dir, 'data')
    images_dir = join(base_dir, 'images')

    def timeit(f, params):
        start_time  = time.time()
        f(*params)
        print(time.time() - start_time)

    timeit(monotonization1, (data_dir, images_dir))
    timeit(monotonization2, (data_dir, images_dir))
