from scipy.interpolate import interp1d
import numpy as np
from fenics import UserExpression

def get_initial_conditions():
    x = (0, 2, 10, 12, 22, 24)
    x = np.array(x) * 3600
    y = (70, 110, 110, 30, 30, 70)
    return interp1d(x, y)


class BC_periodic(UserExpression):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ic = get_initial_conditions()
        self.t = 0
        self.rho = args[0]

    def eval(self, value, x):
        value[0] = self.rho * self.ic(self.t)
    
    def update_t(self, t):
        self.t = t if t <= 24 * 3600 else t - 24 * 3600
    
    def value_shape(self):
        return ()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from os.path import join


    f = get_initial_conditions()

    x = np.linspace(0, 24*3600, 25)
    y = f(x)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    plt.plot(x, y)

    plt.xticks(range(0, 25*3600, 4*3600))

    plt.xlabel(r'$t$')
    plt.ylabel(r'$q_{out}$')
    plt.xlim(x[0], x[-1])
    #plt.legend()
    plt.grid()
    plt.savefig(join('images', 'pipes', 'periodic_bc.pdf'), transparent=True)
