from os.path import dirname, join, exists
from os import makedirs
import numpy as np
from multiprocessing import Pool
import sys
import nonlinear_viscosity
import viscosity1
import viscosity2
import weight
sys.path.append('calculations')
import plotting


def store(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest, calculate, fname):
    mesh_size = 200
    tau = 0.0025
    T = 0.6
    ts2store = [0.0, 0.2, 0.4, 0.6]
    try:
        x, u, t, m = calculate(mesh_size, tau, degree, T, ts2store, sigma, alfa, gamma, kappa, ntest)
    except:
        return
    np.savez_compressed(fname, x=x, u=u, t=t, m=m, ts=np.array(ts2store))


def plot(mesh_size, tau, degree, alfa, gamma, kappa, ntest, fname, images_dir):
    if exists(fname):
        npz = np.load(fname)

        plotting.base([npz['x'] for i in range(npz['ts'].size)], npz['u'], r'$x$', r'$u$', [f'$t={i}$' for i in npz['ts']],
                  join(images_dir, f'u_ms{mesh_size}_tau{tau}_d{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png'))
        
        plotting.base([npz['t']], [npz['m']], r'$t$', r'$m$', [],
                      join(images_dir, f'm_ms{mesh_size}_tau{tau}_d{degree}_alfa{alfa}_gamma{gamma}_kappa{kappa}_ntest{ntest}.png'))


def process_nonlinear_viscosity(mesh_size, tau, p, T, ts2store, a, y, test, data_dir, images_dir):
        try:
            x, u, x_e, u_e, t, err = nonlinear_viscosity.calculate(mesh_size, tau, p, T, ts2store, a, y, test)
        except:
            return
        ts2store = np.array(ts2store)
        np.savez_compressed(join(data_dir, f'a{a}_y{y}_p{p}_test{test}.npz'), x=x, u=u, x_e=x_e, u_e=u_e, t=t, err=err, ts=ts2store)

        # plotting.base([x for i in range(ts2store.size)], u, r'$x$', r'$u$', [f'$t={i}$' for i in ts2store],
        #               join(images_dir, f'a{a}_y{y}_p{p}_test{test}.png'))
        
        # plotting.base([t], [m], r'$t$', r'$m$', [],
        #               join(images_dir, f'm_a{a}_y{y}_p{p}_test{test}.png'))


def experiments_nonlinear_viscosity(data_dir, images_dir):
    name = 'nonlinear_viscosity'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)
    
    mesh_size = 200
    tau = 0.0025
    T = 0.6
    ts2store = [0.2, 0.4, 0.6]

    tests = [2]

    degrees = [1, 2] #list(range(1, 4))
    alfas = [0, 1, 128, 512, 1024]#[0, 1, 16, 32] #[0] + [2**i for i in range(8)]
    gammas = list(range(4))

    param_combinations = [
        (mesh_size, tau, p, T, ts2store, a, y, test, data_dir, images_dir)
        for p in degrees
        for a in alfas
        for y in gammas
        for test in tests
        if not (a == 0 and y != 0)
    ]

    with Pool() as p:
        res1 = p.starmap_async(process_nonlinear_viscosity, param_combinations)
        res1.wait()


def process(mesh_size, tau, p, T, ts2store, k, test, data_dir, images_dir, calculate):
        try:
            x, u, x_e, u_e, t, err = calculate(mesh_size, tau, p, T, ts2store, k, test)
        except:
            return
        ts2store = np.array(ts2store)
        np.savez_compressed(join(data_dir, f'k{k}_test{test}.npz'), x=x, u=u, x_e=x_e, u_e=u_e, t=t, err=err, ts=ts2store)

        # plotting.base([x for i in range(ts2store.size)], u, r'$x$', r'$u$', [f'$t={i}$' for i in ts2store],
        #               join(images_dir, f'k{k}_p{p}_test{test}.png'))
        
        # plotting.base([t], [m], r'$t$', r'$m$', [],
        #               join(images_dir, f'm_k{k}_p{p}_test{test}.png'))


def experiments(kappas, calculate, data_dir, images_dir):
    mesh_size = 200
    tau = 0.0025
    T = 0.6
    ts2store = [0.2, 0.4, 0.6]

    #tests = list(range(1, 4))
    tests = [2]

    degrees = [1]
    #kappas = [0.02, 0.04, 0.08] #[0] + [0.01 * 2**i for i in range(8)]

    param_combinations = [
        (mesh_size, tau, p, T, ts2store, k, test, data_dir, images_dir, calculate)
        for p in degrees
        for k in kappas
        for test in tests
    ]

    with Pool() as p:
        res1 = p.starmap_async(process, param_combinations)
        res1.wait()


def test_viscosity1(data_dir, images_dir):
    name = 'viscosity1'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)
    kappas = [0.16, 1.28, 2.56]
    experiments(kappas, viscosity1.calculate, data_dir, images_dir)


def test_viscosity2(data_dir, images_dir):
    name = 'viscosity2'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)
    kappas = [0.02, 0.04, 0.08]
    experiments(kappas, viscosity2.calculate, data_dir, images_dir)


def test_weight(data_dir, images_dir):
    name = 'weight'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)
    makedirs(images_dir, exist_ok=True)
    kappas = [0.16, 1.28, 2.56]
    experiments(kappas, weight.calculate, data_dir, images_dir)

if __name__ == '__main__':
    import time
    base_dir = dirname(__file__)
    data_dir = join('data', 'monotonization', 'nonlinear')
    images_dir = join('images', 'monotonization', 'nonlinear')

    def timeit(f, params):
        start_time  = time.time()
        f(*params)
        print(time.time() - start_time)
    
    # timeit(test_viscosity1, (data_dir, images_dir))
    # timeit(test_viscosity2, (data_dir, images_dir))
    # timeit(test_weight, (data_dir, images_dir))

    timeit(experiments_nonlinear_viscosity, (data_dir, images_dir))
