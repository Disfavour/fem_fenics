import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


dpi=300
lines = ['-', '--', ':', '-.']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']


def dif_times(xs, yss, ylabel, labels, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    lines = ['-', ':', '--', '-.']  # увеличение монотонности
    if len(yss[0]) == 3:
        lines = ['-', '--', ':']

    for ys, c, label in zip(yss, colors, labels):       # t -> color
        plt.plot(xs[0], ys[0], lines[0]+c, label=label) # sigma -> line type
        for x, y, l in zip(xs[1:], ys[1:], lines[1:]):
            plt.plot(x, y, l+c)

    plt.xlabel(r'$x$')
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname, transparent=True)
    plt.close()


def dif_params(ts, yss, ylabel, labels, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    for t, ys, c, label in zip(ts, yss, colors, labels):    # M/tau -> color
        plt.plot(t, ys[0], lines[0]+c, label=label)         # sigma -> line type
        for y, l in zip(ys[1:], lines[1:]):
            plt.plot(t, y, l+c)

    plt.xlabel(r'$t$')
    plt.ylabel(ylabel)
    plt.xlim(ts[0][0], ts[0][-1])
    plt.legend()
    plt.grid()
    plt.savefig(fname, transparent=True)
    plt.close()


# def dif_times(x, ys, x_e, ys_e, ylabel, names, fname):
#     plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
#     for y, y_e, c, lbl in zip(ys, ys_e, colors, names):
#         plt.plot(x, y, f'--{c}')
#         plt.plot(x_e, y_e, f'-{c}', label=lbl)
#     plt.xlabel(r'$x$')
#     plt.ylabel(ylabel)
#     plt.xlim(x[0], x[-1])
#     plt.legend()
#     plt.grid()
#     plt.savefig(fname)
#     plt.close()


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
    plt.savefig(fname, transparent=True)
    plt.close()


def base(xs, ys, xlabel, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for x, y, c in zip(xs, ys, colors):
        plt.plot(x, y, c)
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
