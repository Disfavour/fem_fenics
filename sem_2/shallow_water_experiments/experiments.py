from os.path import dirname, join
from calcs.shallow_water_1d import shallow_water_1d
from plotting import base_plots


base_dir = dirname(__file__)
plots_dir = join(base_dir, 'plots')


def run_test_problem(sigmas, tau, mesh_size, T, time_moments, critical_time, degree):
    ts = None
    Es = []

    for s in sigmas:
        data = shallow_water_1d(s, tau, mesh_size, T, time_moments, critical_time, degree)
    
        x, h_numerical, u_numerical, h_exact, u_exact, t_L2, h_L2, u_L2, time_steps, m, E = data

        ts = time_steps
        Es.append(E)

        # for i, t in enumerate(time_moments):
        #     base_plots.h(x, h_numerical[i], h_exact[i], plots_dir, f'h_t{t}_s{s}_tau{tau}_ms{mesh_size}')
        #     base_plots.u(x, u_numerical[i], u_exact[i], plots_dir, f'u_t{t}_s{s}_tau{tau}_ms{mesh_size}')

        base_plots.h_at_different_moments(x, h_numerical, h_exact, time_moments, plots_dir, f'h_moments_s{s}_tau{tau}_ms{mesh_size}')
        base_plots.u_at_different_moments(x, u_numerical, u_exact, time_moments, plots_dir, f'u_moments_s{s}_tau{tau}_ms{mesh_size}')

        #base_plots.h_L2(t_L2, h_L2, plots_dir, f'h_L2_s{s}_tau{tau}_ms{mesh_size}')
        #base_plots.u_L2(t_L2, u_L2, plots_dir, f'u_L2_s{s}_tau{tau}_ms{mesh_size}')

    #base_plots.E_different_s(ts, Es, sigmas, plots_dir, f'E_tau{tau}_ms{mesh_size}')


if __name__ == '__main__':
    # sigma
    sigmas = [1]#[0.75, 1, 1.5, 2]
    tau = 0.005
    mesh_size = 200
    # Задача решается от t=0 до T
    T = 1
    # моменты времени для которых будут построены графики. Должно быть 3
    time_moments = [0.1, 0.4, 0.7, 1]
    # только до этого критического времени в графиках будет рисоваться аналитическое решение
    critical_time = 0.5 # подобрано
    degree = 1

    # Графики в plots/{pdf и png}

    # Линия из точек - аналитическое решение

    # Графики на момент времени имеют следующие названия {u или h}_t{момент времени}_s{s}_tau{tau}_ms{mesh_size}
    
    # Графики на разные моменты времени имеют следующие названия {h или u}_moments_s{s}_tau{tau}_ms{mesh_size}
    # поддерживается разная длина списка time_moments (будет нарисовано столько моментов, сколько длина time_moments)

    # Графики нормы L2 абсолютной ошибки имеют следующие названия {h или u}_L2_s{s}_tau{tau}_ms{mesh_size}

    # Графики энергии при разных сигма имеют следующие названия E_tau{tau}_ms{mesh_size}
    # поддерживается разная длина списка sigmas (будет нарисовано столько линий, сколько длина sigmas)

    run_test_problem(sigmas, tau, mesh_size, T, time_moments, critical_time, degree)
