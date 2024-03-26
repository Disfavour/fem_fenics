from scipy.integrate import quad
from scipy.constants import g
import numpy as np


hl, hr = 10, 1
h1 = 3.9618
u1 = g * (h1 * h1 - hr * hr) / (2*g*(h1 + hr)*h1*hr) ** 0.5
D1 = - (g * hl) ** 0.5
D2 = u1 - (g*h1) ** 0.5
D3 = u1 * h1 / (h1 - hr)
domain_size = 5


def integrate(x, y):
    return 0.5*(y[:-1]+y[1:]).sum() * (x[1] - x[0])


def calculate_E(x, h, u):
    return integrate(x, 0.5 * (h * u ** 2 + g * h ** 2))


def calculate_m(x, h):
    return integrate(x, h)


def calculate_h(x, t):
    'x[0] < D1*t ? hl : (x[0] < D2*t ? 1/(9*g) * pow(2*sqrt(g*hl) - x[0]/t, 2) : (x[0] < D3*t ? h1 : hr))'
    if x < D1*t:
        return hl
    elif x < D2*t:
        return 1/(9*g) * (2*(g*hl) ** 0.5 - x/t) ** 2
    elif x < D3*t:
        return h1
    else:
        return hr


def calculate_u(x, t):
    'x[0] < D1*t ? 0  : (x[0] < D2*t ? 1./3 * (2*sqrt(g*hl) + 2*x[0]/t)        : (x[0] < D3*t ? u1 : 0 ))'
    if x < D1*t:
        return 0
    elif x < D2*t:
        return 1/3 * (2*(g*hl)**0.5 + 2*x/t)
    elif x < D3*t:
        return u1
    else:
        return 0


def calculate_m_quad(x, t):
    return calculate_h(x, t)


def calculate_E_quad(x, t):
    h = calculate_h(x, t)
    u = calculate_u(x, t)
    return 0.5 * (h * u ** 2 + g * h ** 2)


def get_solution(xs, ts, time_moments):
    hs, us, ms, Es, ms_quad, Es_quad = [], [], [], [], [], []
    for t in ts:
        h, u = [], []
        for x in xs:
            h.append(calculate_h(x, t))
            u.append(calculate_u(x, t))
        h, u = np.array(h), np.array(u)
        m, E = calculate_m(xs, h), calculate_E(xs, h, u)
        m_quad = quad(calculate_m_quad, -domain_size, domain_size, (t,), limit=200)[0]
        E_quad = quad(calculate_E_quad, -domain_size, domain_size, (t,), limit=200)[0]

        if np.isclose(t, time_moments).any():
            hs.append(h)
            us.append(u)

        ms.append(m)
        Es.append(E)
        ms_quad.append(m_quad)
        Es_quad.append(E_quad)

    return hs, us, ms, Es, ms_quad, Es_quad


if __name__ == '__main__':
    pass
