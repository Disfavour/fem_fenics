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


def store_data(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, calculate, fname):
    try:
        x, u = calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
    except:
        return
    np.save(fname, np.vstack((x, u)))



def calculate_and_plot(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, calculate, fname):
    try:
        x, u = calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
    except:
        return
    plotting.base((x, x, x), u, r'$x$', r'$u$', [f'$t={i}$' for i in ts2store], fname)
    

def experiments(mesh_sizes, taus, degrees, T, sigma, alfas, gammas, kappas, ntests, calculate, data_dir, images_dir):
    ts = np.arange(0, T+taus[-1]/2, taus[-1])
    #ts2store = [ts[ts.shape[0]//3], ts[2*ts.shape[0]//3], ts[-1]]
    ts2store = [ts[0], ts[ts.shape[0]//2], ts[-1]]

    param_combinations = [
        (mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
        for mesh_size in mesh_sizes
        for tau in taus
        for degree in degrees
        for alfa in alfas
        for gamma in gammas
        for kappa in kappas
        for ntest in ntests
        if not exists(join(images_dir, f'ms{mesh_size}_tau{tau}_degree{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png')) 
        and not (alfa == 0 and gamma != 0)
    ]

    fnames = [
        join(images_dir, f'ms{mesh_size}_tau{tau}_degree{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png')
        for (mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest) in param_combinations
    ]

    params = [
        (*params, calculate, fname) for params, fname in zip(param_combinations, fnames)
    ]

    # for p in params:
    #     calculate_and_plot(*p)
    
    with Pool() as p:
        res1 = p.starmap_async(calculate_and_plot, params)
        res1.wait()


def monotonization_degree(data_dir, images_dir):
    name = 'monotonization_degree'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    degrees = list(range(1, 4))
    T = 7.0
    sigma = 0.5
    alfas = [0]
    gammas = [0]
    kappas = [0.5]
    ntests = (1, 2)

    experiments(mesh_sizes, taus, degrees, T, sigma, alfas, gammas, kappas, ntests, monotonization_all.calculate, data_dir, images_dir)


def monotonization_artificial_viscosity(data_dir, images_dir):
    name = 'monotonization_artificial_viscosity'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    degrees = list(range(1, 4))
    T = 7.0
    sigma = 0.5
    alfas = [0] + [2**i for i in range(8)]
    gammas = list(range(4))
    kappas = [0.5]
    ntests = (1, 2)

    experiments(mesh_sizes, taus, degrees, T, sigma, alfas, gammas, kappas, ntests, monotonization_all.calculate, data_dir, images_dir)


def monotonization_kappa(data_dir, images_dir):
    name = 'monotonization_kappa'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_sizes = (800, 400, 200)
    taus = (0.005, 0.01, 0.02)
    degrees = list(range(1, 4))
    T = 7.0
    sigma = 0.5
    alfas = [0]
    gammas = [0]
    kappas = [0.5] + [round(0.5 + 0.01 * 2**i, 2) for i in range(8)]
    ntests = (1, 2)

    experiments(mesh_sizes, taus, degrees, T, sigma, alfas, gammas, kappas, ntests, monotonization_all.calculate, data_dir, images_dir)


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
    alfas = [0] + [2**i for i in range(8)]
    gammas = list(range(4))
    kappas = [0.5] + [round(0.5 + 0.01 * 2**i, 2) for i in range(8)]
    ntests = (1, 2)

    experiments(mesh_sizes, taus, degrees, T, sigma, alfas, gammas, kappas, ntests, monotonization_all.calculate, data_dir, images_dir)


if __name__ == '__main__':
    import time
    base_dir = dirname(__file__)
    data_dir = join(base_dir, 'data')
    images_dir = join(base_dir, 'images')

    def timeit(f, params):
        start_time  = time.time()
        f(*params)
        print(time.time() - start_time)

    timeit(monotonization_degree, (data_dir, images_dir))
    timeit(monotonization_artificial_viscosity, (data_dir, images_dir))
    timeit(monotonization_kappa, (data_dir, images_dir))
        
    #timeit(monotonization_full, (data_dir, images_dir))
