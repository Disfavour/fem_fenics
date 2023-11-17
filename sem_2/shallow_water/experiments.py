from os.path import dirname, join
import matplotlib.pyplot as plt
import numpy as np
from calcs.shallow_water_1d import shallow_water_1d
from plotting.h import h
from plotting.u import u

base_dir = dirname(__file__)
plots = join(base_dir, 'plots')
data = join(base_dir, 'data')


def run_test_problem(s, tau, mesh_size, T, time_to_store_data, time_critical, degree):
    nt = len(time_to_store_data)

    n_with_exact = 0
    for i in time_to_store_data:
        if i < time_critical:
            n_with_exact += 1

    res = shallow_water_1d(s, tau, mesh_size, T, time_to_store_data, time_critical, degree)

    #fname = f's{s}_tau{tau}_ms{mesh_size}'
    #np.save(join(data, fname), res)

    for j, t in enumerate(time_to_store_data):
        i = j + 1
        if nt + i + nt + n_with_exact < res.shape[0]:
            h(res[0], res[i], res[i + 2*nt], plots, f'h_t{t}_s{s}_tau{tau}_ms{mesh_size}')
            u(res[0], res[nt + i], res[nt + i + nt + n_with_exact], plots, f'u_t{t}_s{s}_tau{tau}_ms{mesh_size}')
        else:
            h(res[0], res[i], None, plots, f'h_t{t}_s{s}_tau{tau}_ms{mesh_size}')
            u(res[0], res[nt + i], None, plots, f'u_t{t}_s{s}_tau{tau}_ms{mesh_size}')


if __name__ == '__main__':
    # sigma
    s = 1
    tau = 0.001
    mesh_size = 1000
    # Задача решается от t=0 до T
    T = 1
    # моменты времени для которых будут построены графики
    time_to_store_data = [0, 0.01, 0.1, 0.3, 0.4, 0.45, 0.7]
    time_to_store_data.sort()
    # только до этого критического времени в графиках будет рисоваться аналитическое решение
    time_critical = 0.5
    degree = 1

    run_test_problem(s, tau, mesh_size, T, time_to_store_data, time_critical, degree)
