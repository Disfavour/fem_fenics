import matplotlib.pyplot as plt
import numpy as np
import os.path

base_dir = os.path.dirname(os.path.dirname(__file__))
plots = os.path.join(base_dir, 'plots')
data = os.path.join(base_dir, 'data')


def special(fnames):
    arrs = []
    for fname in fnames:
        arrs.append(np.load(fname))

    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

    plt.plot(arrs[0][0], arrs[0][2], color='blue', linestyle='solid', label=r'$\tau=0.01$')

    plt.plot(arrs[1][0], arrs[1][2], color='red', linestyle="solid", label=r"$\tau=0.005$")

    plt.plot(arrs[2][0], arrs[2][2], color='green', linestyle="solid", label=r"$\tau=0.0025$")

    plt.plot(arrs[3][0], arrs[3][2], color='blue', linestyle='solid', label=r'$\tau=0.01$')

    plt.plot(arrs[4][0], arrs[4][2], color='red', linestyle="solid", label=r"$\tau=0.005$")

    plt.plot(arrs[5][0], arrs[5][2], color='green', linestyle="solid", label=r"$\tau=0.0025$")

    plt.legend()
    plt.grid()
    #plt.xlim(arr_t1[0][0], arr_t1[0][-1])
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$\varrho$')

    for i in ['pdf', 'png']:
        plt.savefig(os.path.join(plots, i, 'special1_' + '.' + i))


if __name__ == '__main__':
    fnames = []
    for prog in ('sym', 'nsym'):
        for tau in ('0.01', '0.005', '0.0025'):
            fnames.append(os.path.join(data, f'{prog}_tau{tau}.npy'))
    special(fnames)
