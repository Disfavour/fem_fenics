import matplotlib.pyplot as plt
import numpy as np
import os.path

base_dir = os.path.dirname(os.path.dirname(__file__))
plots = os.path.join(base_dir, 'plots')
data = os.path.join(base_dir, 'data')


def fig4(fnamet1, fnamet2, fnamet3, resname):
    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    lws = [1.5, 1.2, 0.9]

    if fnamet1 is not None:
        arr_t1 = np.load(fnamet1)
        plt.plot(arr_t1[10], arr_t1[11], ':b', lw=lws[0])
        plt.plot(arr_t1[10], arr_t1[12], ':r', lw=lws[0])
        plt.plot(arr_t1[10], arr_t1[13], ':g', lw=lws[0])

    if fnamet2 is not None:
        arr_t2 = np.load(fnamet2)
        plt.plot(arr_t2[10], arr_t2[11], '--b', lw=lws[1])
        plt.plot(arr_t2[10], arr_t2[12], '--r', lw=lws[1])
        plt.plot(arr_t2[10], arr_t2[13], '--g', lw=lws[1])

    if fnamet3 is not None:
        arr_t3 = np.load(fnamet3)
        plt.plot(arr_t3[10], arr_t3[11], '-b', lw=lws[2], label=r"$t=1$")
        plt.plot(arr_t3[10], arr_t3[12], '-r', lw=lws[2], label=r"$t=2$")
        plt.plot(arr_t3[10], arr_t3[13], '-g', lw=lws[2], label=r"$t=3$")
        plt.xlim(arr_t3[10][0], arr_t3[10][-1])

    plt.legend()
    plt.grid()
    #plt.xlim(arr_t3[10][0], arr_t3[10][-1])
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$\varrho$')

    for i in ['pdf', 'png']:
        plt.savefig(os.path.join(plots, i, 'fig4_' + resname + '.' + i), transparent=True)


if __name__ == '__main__':
    # for prog in ('sym', 'nsym'):
    #     for ms in (100, 200):
    #         fnames = []
    #         for tau in ('0.01', '0.005', '0.0025'):
    #             fnames.append(os.path.join(data, f'{prog}_tau{tau}_ms{ms}.npy'))
    #         fig4(*fnames, f'{prog}_ms{ms}')

    for ms in (100, 200):
        for s in [0, 0.25, 0.5, 0.75, 1, 1.5, 2]:
            fnames = []
            for tau in ('0.01', '0.005', '0.0025'):
                f = os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy')
                f = f if os.path.isfile(f) else None
                fnames.append(f)
            if all(fnames):
                fig4(*fnames, f'w_s{s}_ms{ms}')

    # w_s1_tau0.00125_ms400

    s = 1
    ms = 400

    fnames = []
    for tau in ('0.0025', '0.00125', '0.000625'):
        f = os.path.join(data, f'w_s{s}_tau{tau}_ms{ms}.npy')
        f = f if os.path.isfile(f) else None
        fnames.append(f)
    if all(fnames):
        fig4(*fnames, f'small_tau_w_s{s}_ms{ms}')
