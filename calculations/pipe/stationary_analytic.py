import numpy as np
from  scipy.optimize import newton


def get_exact(xs, D, A, Rs, T, f, m, P_left):
    gamma = Rs * T * m**2 / A**2
    B = f * m**2 * Rs * T / (2 * D * A**2)
    C = P_left**2 / 2 - gamma * np.log(P_left)

    f = lambda p, x: p**2/2 - gamma * np.log(p) + B*x - C
    fprime = lambda p, x: p - gamma / p
    fprime2 = lambda p, x: 1 + gamma / p**2

    exact = [newton(f, P_left, fprime, (xs[0],), fprime2=fprime2)]

    for x in xs[1:]:
        exact.append(newton(f, exact[-1], fprime, (x,),  fprime2=fprime2))

    return np.array(exact)
