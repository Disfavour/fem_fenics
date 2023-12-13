from fenics import *
import math
import matplotlib.pyplot as plt
import numpy as np

set_log_level(LogLevel.WARNING)

a = 1
b = 0.5
c = 2
pi = math.pi
alpha = 1e-3#1e-1
gamma = 1e-6#alpha/2
resW = 1000

degree = 1

mesh = UnitIntervalMesh(resW)

F = FunctionSpace(mesh, 'P', degree)

source = '- cos(pi * x[0]) * (a*pi + (c/pi)) + ' \
        'sin(pi * x[0]) * b + ' \
        'x[0] * (2*a + b) +' \
        '(pow(x[0],2)) * ((c/2) - b)-' \
        '(pow(x[0], 3)) * c/3 -' \
        'a - c * ((1/pi) + (1/6))  '
source_f_ = Expression(source, degree=degree, a = a, b = b, c = c, pi = pi)
source_f = project(source_f_, F)

phi_obs = '- (1/pi)*cos(pi * x[0]) + (pow(x[0], 2))/2 - (pow(x[0], 3))/3 - ((1/pi) + (1/6))'
exact = Expression(phi_obs, degree=degree, pi = pi, eps=1.0e-9)
u_ex = project(exact, F)

manage = Function(F)
manage_next = Function(F)

# 1
phi = TrialFunction(F)
phi_t = TestFunction(F)

F1 = (a * inner(grad(phi), grad(phi_t)) * dx) \
    + (b * Dx(phi,0) * phi_t * dx) \
    + (c * dot(phi, phi_t) * dx) \
    - (dot(source_f + manage, phi_t) * dx)

a1 = lhs(F1)
L1 = rhs(F1)

phi = Function(F)

# 2
phi_s = TrialFunction(F)
phi_s_t = TestFunction(F)

F2 = (a * inner(grad(phi_s), grad(phi_s_t)) * dx) \
        - (b * Dx(phi_s,0) * phi_s_t * dx) \
        + (c * phi_s * phi_s_t * dx) \
        - ((phi - u_ex)* phi_s_t * dx)


a2 = lhs(F2)
L2 = rhs(F2)

phi_s = Function(F)

manage.assign(project(Constant(0), F))

for iter in range(10):
    solve(a1 == L1, phi)
    solve(a2 == L2, phi_s)

    print('phi, u_ex', errornorm(phi, u_ex, norm_type='L2', degree_rise=0))

    manage_next.assign(manage + gamma * (alpha * manage + phi_s))

    print('manage', errornorm(manage_next, manage, norm_type='L2', degree_rise=0))

    manage.assign(manage_next)

    k = alpha * np.sqrt((manage.compute_vertex_values()**2).sum()) + np.sqrt(((phi.compute_vertex_values() - u_ex.compute_vertex_values())**2).sum())
    print('k=', k)
    print()


#     plt.figure()
#     plot(phi)
# plt.show()


# q1 = project(Constant(1), F)
# q2 = project(Constant(1.333), F)

# q3 = Function(F)
# q2.assign(q1 + gamma * (alpha * q1 + q2))

# plot(q2)
# plt.show()
