import matplotlib.pyplot as plt
import numpy as np
import os.path
import matplotlib.ticker as mtick

base_dir = os.path.dirname(os.path.dirname(__file__))
plots = os.path.join(base_dir, 'plots')
data = os.path.join(base_dir, 'data')


def fig_dif_sigma(fnames, resname):
    arrs = []

    for fname in fnames:
        arrs.append(np.load(fname))

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    lws = [0.9, 1.2, 1.5]

    plt.plot(arrs[0][0], arrs[0][2], '-b', lw=lws[0], label=r'$\sigma=0.5$')

    plt.plot(arrs[1][0], arrs[1][2], '-r', lw=lws[0], label=r'$\sigma=0.75$')
    plt.plot(arrs[2][0], arrs[2][2], '-g', lw=lws[0], label=r'$\sigma=1$')

    plt.plot(arrs[3][0], arrs[3][2], '-c', lw=lws[0], label=r'$\sigma=1.5$')
    plt.plot(arrs[4][0], arrs[4][2], '-m', lw=lws[0], label=r'$\sigma=2$')

    # if fnamet1 is not None:
    #     arr_t1 = np.load(fnamet1)
    #     plt.plot(arr_t1[0], arr_t1[2], ':b', lw=1.5, label=r'$\sigma=0.5$')

    # if fnamet2 is not None:
    #     arr_t2 = np.load(fnamet2)
    #     plt.plot(arr_t2[0], arr_t2[2], '--r', lw=1.2, label=r"$\sigma=0.75$")

    # if fnamet3 is not None:
    #     arr_t3 = np.load(fnamet3)
    #     plt.plot(arr_t3[0], arr_t3[2], '-g', lw=0.9, label=r"$\sigma=1$")
    #     plt.xlim(arr_t3[0][0], arr_t3[0][-1])

    #plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))

    plt.legend()
    plt.grid()
    plt.xlim(arrs[0][0][0], arrs[0][0][-1])
    plt.xlabel(r"$t$")
    plt.ylabel(r"$E$")

    for i in ['pdf', 'png']:
        plt.savefig(os.path.join(plots, i, 'fig_dif_sigma_' + resname + '.' + i), transparent=True)


if __name__ == '__main__':
    # for prog in ('sym', 'nsym'):
    #     for ms in (100, 200):
    #         fnames = []
    #         for tau in ('0.01', '0.005', '0.0025'):
    #             fnames.append(os.path.join(data, f'{prog}_tau{tau}_ms{ms}.npy'))
    #         fig11(*fnames, f'{prog}_ms{ms}')

    for ms in (100, 200):
        for tau in ('0.01', '0.005', '0.0025'):
            fnames = []
            for s in [0.5, 0.75, 1, 1.5, 2]:
                f = os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy')
                f = f if os.path.isfile(f) else None
                fnames.append(f)
            if all(fnames):
                fig_dif_sigma(fnames, f'w_t{tau}_ms{ms}')
