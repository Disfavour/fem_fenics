from multiprocessing import Pool
from symmetric import symmetric
from nonsymmetric import nonsymmetric
from weighted import weighted


if __name__ == '__main__':
    mesh_sizes = [200, 100]
    taus = [0.0025, 0.005, 0.01]
    sigmas = [1.5, 2]
    params = ((s, tau, ms, None, True) for ms in mesh_sizes for tau in taus for s in sigmas)
    # print(*params)
    # exit()

    with Pool() as p:
        res = p.starmap_async(weighted, params)
        res.wait()
