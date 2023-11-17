import matplotlib.pyplot as plt

dpi = 300


def density_at_different_moments(x, r1, r2, r3, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(x, r1, label=r'$t=1$')
    plt.plot(x, r2, label=r'$t=2$')
    plt.plot(x, r3, label=r'$t=3$')

    plt.xlim(x[0], x[-1])
    plt.xlabel(r'$x_1$')
    plt.ylabel(r'$\varrho$')
    plt.legend()
    plt.grid()

    plt.savefig(filename)
    plt.close()


def E(t, E, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(t, E)

    plt.xlim(t[0], t[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.grid()

    plt.savefig(filename)
    plt.close()


def E_compare(t, E1, E2, E3, E4, label1, label2, label3, label4, filename):
    plt.figure(figsize=(6.4, 3.6), dpi=dpi, tight_layout=True)

    plt.plot(t, E1, label=label1)
    plt.plot(t, E2, label=label2)
    plt.plot(t, E3, label=label3)
    plt.plot(t, E4, label=label4)

    plt.xlim(t[0], t[-1])
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    plt.grid()
    plt.legend()

    plt.savefig(filename)
    plt.close()


if __name__ == '__main__':
    pass
