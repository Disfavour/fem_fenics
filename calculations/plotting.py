import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


dpi=300
lines = ['-', '--', '-.', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']


def dif_times(x, ys, x_e, ys_e, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for y, y_e, c, lbl in zip(ys, ys_e, colors, names):
        plt.plot(x, y, f'--{c}')
        plt.plot(x_e, y_e, f'-{c}', label=lbl)
    plt.xlabel(r'$x$')
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()


def dif_times_dif_param(x, yss, x_e, y_e, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for ys, l in zip(yss, lines[1:]):
        for y, c in zip(ys, colors):
            plt.plot(x, y, l+c)
    
    for y, c, n in zip(y_e, colors, names):
        plt.plot(x_e, y, lines[0]+c, label=n)

    plt.xlabel(r'$x$')
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()


def dif_times_dif_param_linear(x, yss, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for ys, c, name in zip(yss, colors, names):
        for i, y in enumerate(ys):
            if i == 0:
                plt.plot(x, y, c, label=name)
            else:
                plt.plot(x, ys.T, c)

    plt.xlabel(r'$x$')
    plt.ylabel(r'$u$')
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname)
    plt.close()


def err_dif_M(t, yss, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for ys, l in zip(yss, lines):
        for y, c in zip(ys, colors):
            plt.plot(t, y, l+c)

    plt.xlabel(r'$t$')
    plt.ylabel(ylabel)
    plt.xlim(t[0], t[-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()


def err_dif_tau(ts, yss, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for ys, l in zip(yss, lines):
        for t, y, c in zip(ts, ys, colors):
            plt.plot(t, y, l+c)

    plt.xlabel(r'$t$')
    plt.ylabel(ylabel)
    plt.xlim(t[0], t[-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()


def base(xs, ys, xlabel, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for x, y in zip(xs, ys):
        plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xs[0][0], xs[0][-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()


def base2(xs, ys, xs2, ys2, xlabel, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for x, y, c in zip(xs, ys, colors):
        plt.plot(x, y, f'-{c}')
    for x, y, c in zip(xs2, ys2, colors):
        plt.plot(x, y, f':{c}')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xs[0][0], xs[0][-1])
    plt.legend(names)
    plt.grid()
    plt.savefig(fname)
    plt.close()


if __name__ == '__main__':
    pass
