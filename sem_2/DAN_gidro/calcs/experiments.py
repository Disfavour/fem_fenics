from multiprocessing import Pool
from symmetric import symmetric
from nonsymmetric import nonsymmetric
from weighted import weighted


if __name__ == '__main__':
    taus = [0.01, 0.005, 0.0025]
    ms = 100
    #params = ((tau, ms, None, True) for tau in taus)
    params = ((s, tau, ms, None, True) for s in [0, 0.25, 0.5, 0.75, 1] for tau in taus)

    with Pool() as p:
        res = p.starmap_async(weighted, params)
        res.wait()
