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


M_s = [100, 200, 400] # [200]
tau_s = [0.005, 0.0025, 0.00125] # [0.0025]
ts2store = [0.2, 0.4, 0.6]

k_s = [0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]


def store(args):
    params = args[:-2]
    fname = args[-2]
    calculate = args[-1]
    try:
        x, u, x_e, u_e, t, err = calculate(*params)
    except:
        return
    np.savez_compressed(fname, x=x, u=u, x_e=x_e, u_e=u_e, t=t, err=err, ts=np.array(params[-1]))


def run_experiments(params):
    
    with Pool() as p:
        res1 = p.starmap_async(store, params)
        res1.wait()
    
    # for p in params:
    #     store(p)


def test_weight(data_dir, images_dir):
    name = 'weight'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)

    params = [
        [[M, tau, k, ts2store, join(data_dir, f'M{M}_tau{tau}_k{k}.npz'), weight.calculate]]
        for M in M_s
        for tau in tau_s
        for k in k_s
    ]

    run_experiments(params)


def test_viscosity1(data_dir, images_dir):
    name = 'viscosity1'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)

    params = [
        [[M, tau, k, ts2store, join(data_dir, f'M{M}_tau{tau}_k{k}.npz'), viscosity1.calculate]]
        for M in M_s
        for tau in tau_s
        for k in k_s
    ]

    run_experiments(params)


def test_viscosity2(data_dir, images_dir):
    name = 'viscosity2'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)

    params = [
        [[M, tau, k, ts2store, join(data_dir, f'M{M}_tau{tau}_k{k}.npz'), viscosity2.calculate]]
        for M in M_s
        for tau in tau_s
        for k in [0.02, 0.04, 0.08]
    ]

    run_experiments(params)


def test_nonlinear_viscosity(data_dir, images_dir):
    name = 'nonlinear_viscosity'
    data_dir = join(data_dir, name)
    images_dir = join(images_dir, name)
    makedirs(data_dir, exist_ok=True)

    a_s = [0, 50, 100, 200, 500, 1000]
    y_s = [0]

    params = [
        [[M, tau, a, y, ts2store, join(data_dir, f'M{M}_tau{tau}_a{a}_y{y}.npz'), nonlinear_viscosity.calculate]]
        for M in M_s
        for tau in tau_s
        for a in a_s
        for y in y_s
    ]

    run_experiments(params)

    a_s = [1]
    y_s = [0, 1, 2, 3]

    params = [
        [[M, tau, a, y, ts2store, join(data_dir, f'M{M}_tau{tau}_a{a}_y{y}.npz'), nonlinear_viscosity.calculate]]
        for M in M_s
        for tau in tau_s
        for a in a_s
        for y in y_s
    ]

    run_experiments(params)

if __name__ == '__main__':
    import time
    base_dir = dirname(__file__)
    data_dir = join('data', 'monotonization', 'nonlinear')
    images_dir = join('images', 'monotonization', 'nonlinear')

    def timeit(f, params):
        start_time  = time.time()
        f(*params)
        print(time.time() - start_time)
    
    timeit(test_weight, (data_dir, images_dir))
    #timeit(test_viscosity1, (data_dir, images_dir))
    #timeit(test_viscosity2, (data_dir, images_dir))
    #timeit(test_nonlinear_viscosity, (data_dir, images_dir))
