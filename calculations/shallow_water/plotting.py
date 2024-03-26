import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


dpi=300
lines = ['-', '--', '-.', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']


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

def dif_times(x, ys, x2, ys2, ylabel, names, fname):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)
    for y, y2, c, lbl in zip(ys, ys2, colors, names):
        plt.plot(x, y, f'-{c}', label=lbl)
        plt.plot(x2, y2, f':{c}')
    plt.xlabel(r'$x$')
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.legend()
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
