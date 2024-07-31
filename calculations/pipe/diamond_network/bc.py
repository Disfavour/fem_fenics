from scipy.interpolate import interp1d
import numpy as np
from fenics import UserExpression

def get_bc():
    x = (0, 1, 3, 6, 8, 12)
    x = np.array(x) * 3600
    y = (100, 100, 200, 200, 80, 80)
    return interp1d(x, y)


class BC_diamond(UserExpression):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bc = get_bc()
        self.t = 0

    def eval(self, value, x):
        value[0] = self.bc(self.t)
    
    def update_t(self, t):
        self.t = t
    
    def value_shape(self):
        return ()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from os.path import join


    f = get_bc()

    x = np.linspace(0, 12*3600, 13)
    y = f(x)

    plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

    plt.plot(x, y)

    plt.xticks(range(0, 13*3600, 2*3600))

    plt.xlabel(r'$t$')
    plt.ylabel(r'$m$')
    plt.xlim(x[0], x[-1])
    #plt.legend()
    plt.grid()
    plt.savefig(join('images', 'pipes', 'diamond_bc.pdf'), transparent=True)
