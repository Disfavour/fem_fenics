from scipy.interpolate import interp1d
import numpy as np
from fenics import UserExpression

def get_bcs():
    x = (0, 4, 12, 20, 24)
    x = np.array(x) * 3600
    y2 = np.array((20, 30, 10, 30, 20))
    y3 = y2 + 20
    return interp1d(x, y2), interp1d(x, y3)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from os.path import join

    bc2, bc3 = get_bcs()

    t = np.linspace(0, 24*3600, 25)
    m2 = bc2(t)
    m3 = bc3(t)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    plt.plot(t, m2, label='2')
    plt.plot(t, m3, label='3')

    plt.xticks(range(0, 25*3600, 4*3600))

    plt.xlabel(r'$t$')
    plt.ylabel(r'$q$')
    plt.xlim(t[0], t[-1])
    plt.legend()
    plt.grid()
    plt.savefig(join('images', 'pipes', 'simple_bc.pdf'), transparent=True)