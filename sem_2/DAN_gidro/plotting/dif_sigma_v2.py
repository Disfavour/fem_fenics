import matplotlib.pyplot as plt
import numpy as np
import os.path
import matplotlib.ticker as mtick

base_dir = os.path.dirname(os.path.dirname(__file__))
plots = os.path.join(base_dir, 'plots')
data = os.path.join(base_dir, 'data')


def fig_dif_sigma(fnames, resname):
    a = []

    for i in fnames:
        a.append(np.load(i))

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    lws = [1.5, 1.0, 0.5]
    lws = [1.5, 1.2, 0.9]

    plt.plot(a[2][0], a[2][2], '-r', lw=lws[2], label=r"$\tau=0.01$")
    plt.plot(a[5][0], a[5][2], '-g', lw=lws[2], label=r"$\tau=0.005$")
    plt.plot(a[8][0], a[8][2], '-b', lw=lws[2], label=r"$\tau=0.0025$")
    plt.plot(a[1][0], a[1][2], '--r', lw=lws[1])
    plt.plot(a[4][0], a[4][2], '--g', lw=lws[1])
    plt.plot(a[7][0], a[7][2], '--b', lw=lws[1])
    plt.plot(a[0][0], a[0][2], ':r', lw=lws[0])
    plt.plot(a[3][0], a[3][2], ':g', lw=lws[0])
    plt.plot(a[6][0], a[6][2], ':b', lw=lws[0])

    # plt.plot(a[0][0], a[0][2], ':r', lw=lws[0])
    # plt.plot(a[1][0], a[1][2], '--r', lw=lws[1])
    # plt.plot(a[2][0], a[2][2], '-r', lw=lws[2], label=r"$\tau=0.01$")
    #
    #
    # plt.plot(a[3][0], a[3][2], ':g', lw=lws[0])
    # plt.plot(a[4][0], a[4][2], '--g', lw=lws[1])
    # plt.plot(a[5][0], a[5][2], '-g', lw=lws[2], label=r"$\tau=0.005$")
    #
    #
    # plt.plot(a[6][0], a[6][2], ':b', lw=lws[0])
    # plt.plot(a[7][0], a[7][2], '--b', lw=lws[1])
    # plt.plot(a[8][0], a[8][2], '-b', lw=lws[2], label=r"$\tau=0.0025$")

    #plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))

    plt.legend()
    plt.grid()
    plt.xlim(a[0][0][0], a[0][0][-1])
    plt.xlabel(r"$t$")
    plt.ylabel(r"$E$")

    for i in ['pdf', 'png']:
        plt.savefig(os.path.join(plots, i, 'fig_dif_sigma_v2_' + resname + '.' + i), transparent=True)


if __name__ == '__main__':
    # for prog in ('sym', 'nsym'):
    #     for ms in (100, 200):
    #         fnames = []
    #         for tau in ('0.01', '0.005', '0.0025'):
    #             fnames.append(os.path.join(data, f'{prog}_tau{tau}_ms{ms}.npy'))
    #         fig11(*fnames, f'{prog}_ms{ms}')

    for ms in (100, 200, 400):
        fnames = []
        for tau in ('0.01', '0.005', '0.0025'):
            for s in [0.5, 0.75, 1]:
                f = os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy')
                f = f if os.path.isfile(f) else None
                fnames.append(f)
        if all(fnames):
            fig_dif_sigma(fnames, f'w_ms{ms}')
