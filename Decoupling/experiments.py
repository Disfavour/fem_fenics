from multiprocessing import Pool
from nonlinear import nonlinear
from jacobi import jacobi
from seidel import seidel


if __name__ == '__main__':
    taus = [0.01, 0.005, 0.0025]
    params_nonlinear = [(1, 100, tau) for tau in taus]
    params = [(k, 1, 100, tau) for k in [1, 2, 5] for tau in taus]

    with Pool() as p:
        res1 = p.starmap_async(nonlinear, params_nonlinear)
        res2 = p.starmap_async(jacobi, params)
        res3 = p.starmap_async(seidel, params)
        res1.wait()
        res2.wait()
        res3.wait()

        #map()
        #print(p.starmap(f, [(1, 1), (2, 2), (3, 3)]))
        # res = p.map_async(g, [1, 2, 3])
        # res1 = p.map_async(g, [4, 5, 6])
        # print(res.get())
        # print(res1.get())

        # res = [p.apply_async(f, t) for t in [(1, 1), (2, 2), (3, 3)]]
        # for i in res:
        #     print(i.get())

