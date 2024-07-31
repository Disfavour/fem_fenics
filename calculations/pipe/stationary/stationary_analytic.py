import numpy as np
from  scipy.optimize import newton


def get_exact(xs, D, A, Rs, T, f, m, P_in):     # m = m_out
    C1 = Rs*T*m**2/A**2
    C2 = f*Rs*T*m*abs(m)/(2*D*A**2)

    eq = lambda p, x: p**2/2 - C1*np.log(abs(p)) + C2*x - P_in**2/2 + C1*np.log(abs(P_in))
    eq_derivative = lambda p, x: p - C1 / p + x - x
    eq_derivative2 = lambda p, x: 1 + C1 / p**2 + x - x

    return newton(eq, np.full(xs.shape, P_in), eq_derivative, (xs,), fprime2=eq_derivative2)
