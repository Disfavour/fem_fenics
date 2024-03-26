from multiprocessing import Pool
from nonlinear import nonlinear
from jacobi import jacobi
from seidel import seidel
from nonlinear_RT import nonlinear_RT
from decoupled_RT import decoupled_RT


if __name__ == '__main__':
    # nonlinear(degree=1, mesh_size=100, tau=0.005, T=5, t_fig4=(1, 2, 3, 4, 5))
    # (1, 100, 0.005)
    # jacobi(K=1, degree=1, mesh_size=100, tau=0.005, T=5, t_fig4=(1, 2, 3, 4, 5))
    # (1, 1, 100, 0.005)
    taus = [0.01, 0.005, 0.0025]
    params_nonlinear = [(1, 100, tau) for tau in taus]
    params = [(k, 3, 100, tau) for k in [1, 2, 5] for tau in taus]

    with Pool() as p:
        res1 = p.starmap_async(nonlinear_RT, params_nonlinear)
        #res2 = p.starmap_async(decoupled_RT, params)
        # res3 = p.starmap_async(seidel, params)
        res1.wait()
        #res2.wait()
        # res3.wait()

        #map()
        #print(p.starmap(f, [(1, 1), (2, 2), (3, 3)]))
        # res = p.map_async(g, [1, 2, 3])
        # res1 = p.map_async(g, [4, 5, 6])
        # print(res.get())
        # print(res1.get())

        # res = [p.apply_async(f, t) for t in [(1, 1), (2, 2), (3, 3)]]
        # for i in res:
        #     print(i.get())

