import matplotlib.pyplot as plt
import numpy as np
import os.path

base_dir = os.path.dirname(os.path.dirname(__file__))
plots = os.path.join(base_dir, 'plots')
data = os.path.join(base_dir, 'data')


def fig2(fname, resname):
    arr = np.load(fname)
    ts, r_00, r_min, r_max = arr[0], arr[7], arr[8], arr[9]

    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

    plt.plot(ts, r_00, '-b', label=r'$\varrho_0$')
    plt.plot(ts, r_min, '-r', label=r'$\varrho_{min}$')
    plt.plot(ts, r_max, '-g', label=r'$\varrho_{max}$')

    plt.legend()
    plt.grid()
    plt.xlim(ts[0], ts[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\varrho$')

    for i in ['pdf', 'png']:
        plt.savefig(os.path.join(plots, i, 'fig2_' + resname + '.' + i))


if __name__ == '__main__':
    for prog in ('sym', 'nsym'):
        for tau in ('0.01', '0.005', '0.0025'):
            for ms in (100, 200):
                fig2(os.path.join(data, f'{prog}_tau{tau}_ms{ms}.npy'), f'{prog}_tau{tau}_ms{ms}')

    ms = 100
    for tau in ('0.01', '0.005', '0.0025'):
        for s in [0, 0.25, 0.5, 0.75, 1]:
            if os.path.isfile(os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy')):
                fig2(os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy'), f'w_s{s}_tau{tau}_ms{ms}')

