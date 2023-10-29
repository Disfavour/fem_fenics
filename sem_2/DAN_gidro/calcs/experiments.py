from multiprocessing import Pool
from symmetric import symmetric
from nonsymmetric import nonsymmetric
from weighted import weighted


if __name__ == '__main__':
    taus = [0.0025, 0.005, 0.01]
    ms = 400
    #params = ((tau, ms, None, True) for tau in taus)
    params = ((s, tau, ms, None, True) for tau in taus for s in [0.5, 0.75, 1])
    # print(*params)

    with Pool() as p:
        res = p.starmap_async(weighted, params)
        res.wait()
