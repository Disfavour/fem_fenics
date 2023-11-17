import matplotlib.pyplot as plt
import numpy as np
from os.path import join


def h(x, h, he, dir, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    lws = [1.5, 1.2]

    plt.plot(x, h, '-b', lw=lws[0], label=r"numerical")

    if he is not None:
        plt.plot(x, he, ':r', lw=lws[1], label=r"exact")

    plt.legend()
    plt.grid()
    plt.xlim(x[0], x[-1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$h$')

    for i in ['pdf', 'png']:
        plt.savefig(join(dir, i, filename), format=i, transparent=True)


if __name__ == '__main__':
    pass
