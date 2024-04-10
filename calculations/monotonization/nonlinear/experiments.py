from os.path import dirname, join, exists
from os import makedirs
import numpy as np
from multiprocessing import Pool
import sys
import nonlinear
sys.path.append('calculations')
import plotting


def store(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, calculate, fname):
    try:
        x, u, xe, ue, t, err = calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
    except:
        return
    np.savez_compressed(fname, x=x, u=u, xe=xe, ue=ue, t=t, err=err, ts=np.array(ts2store))


def plot(mesh_size, tau, degree, alfa, gamma, kappa, ntest, fname, images_dir):
    if exists(fname):
        npz = np.load(fname)

        plotting.dif_times(npz['x'], npz['u'],
                           npz['xe'], npz['ue'],
                           r'$u$',
                           [f'$t={i}$' for i in npz['ts']],
                           join(images_dir, f'u_ms{mesh_size}_tau{tau}_d{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png'))
        
        plotting.base([npz['t']], [npz['err']], r'$t$', r'$\varepsilon$', [],
                      join(images_dir, f'err_ms{mesh_size}_tau{tau}_d{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png'))


def experiments(mesh_size, tau, degrees, T, ts2store, sigma, alfas, gammas, kappas, ntests, calculate, data_dir, images_dir):
    param_combinations = [
        (mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
        for degree in degrees
        for alfa in alfas
        for gamma in gammas
        for kappa in kappas
        for ntest in ntests
        if not (alfa == 0 and gamma != 0)
    ]

    fnames_data = [
        join(data_dir, f'ms{mesh_size}_tau{tau}_d{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.npz')
        for (mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest) in param_combinations
    ]

    params_data = [
        (*params, calculate, fname) for params, fname in zip(param_combinations, fnames_data)
    ]
    
    with Pool() as p:
        res1 = p.starmap_async(store, params_data)
        res1.wait()

    param_combinations_images = [
        (mesh_size, tau, degree, alfa, gamma, kappa, ntest, fname, images_dir)
        for (mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, calculate, fname) in params_data
    ]

    with Pool() as p:
        res1 = p.starmap_async(plot, param_combinations_images)
        res1.wait()


def viscosity(data_dir, images_dir):
    name = 'viscosity'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_size = 200
    tau = 0.0025
    degrees = list(range(1, 4))
    T = 0.6
    ts2store = [0.0, 0.2, 0.4, 0.6]
    sigma = 0.5
    alfas = [0] + [2**i for i in range(8)]
    gammas = list(range(4))
    kappas = [0]
    ntests = (1, 2)

    experiments(mesh_size, tau, degrees, T, ts2store, sigma, alfas, gammas, kappas, ntests, nonlinear.calculate, data_dir, images_dir)


def weight(data_dir, images_dir):
    name = 'weight'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)

    mesh_size = 200
    tau = 0.0025
    degrees = list(range(1, 4))
    T = 0.6
    ts2store = [0.0, 0.2, 0.4, 0.6]
    sigma = 0.5
    alfas = [0]
    gammas = [0]
    kappas = [0] + [0.01 * 2**i for i in range(8)]
    ntests = (1, 2)

    experiments(mesh_size, tau, degrees, T, ts2store, sigma, alfas, gammas, kappas, ntests, nonlinear.calculate, data_dir, images_dir)


if __name__ == '__main__':
    import time
    base_dir = dirname(__file__)
    data_dir = join('data', 'monotonization', 'nonlinear')
    images_dir = join('images', 'monotonization', 'nonlinear')

    def timeit(f, params):
        start_time  = time.time()
        f(*params)
        print(time.time() - start_time)

    timeit(viscosity, (data_dir, images_dir))
    timeit(weight, (data_dir, images_dir))
