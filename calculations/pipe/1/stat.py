from fenics import *
import numpy as np
from scipy.constants import pi, R
import matplotlib.pyplot as plt
from  scipy.optimize import newton


set_log_level(LogLevel.WARNING)

P_in, P_out, m_in, m_out = [], [], [], []

mesh_size = 100

L = 1e5
T = 278

D = 0.6
A = pi * D**2 / 4

eps = 0.000617
Re = 14000
f = (-2*np.log10(eps/D/3.7 - 4.518/Re*np.log10(6.9/Re + (eps/D/3.7)**1.11))) ** -2

S = 0.6
M_air = 28.964917 / 1000
M = S * M_air
Rs = R / M

P_left = 5e6
#rho = P_left / (Rs * T)
rho = 0.73

m = rho * 70

mesh = IntervalMesh(mesh_size, 0, L)

P = FunctionSpace(mesh, 'P', 1)

bc = DirichletBC(P, Constant(P_left), 'on_boundary && near(x[0], 0)')

p = Function(P)
pt = TestFunction(P)

F = (p.dx(0) + Rs*T*m**2/A**2 * (p**-1).dx(0) + Rs*T*f*m**2/(2*D*A**2*p)) * pt*dx

p.assign(project(Constant(P_left), P))

solve(F == 0, p, bc)

gamma = Rs * T * m**2 / A**2
B = f * m**2 * Rs * T / (2 * D * A**2)
C = P_left**2 / 2 - gamma * np.log(P_left)
print(gamma, B, C)

xs = mesh.coordinates().flatten()

res = []
for x in xs:
    res.append(newton(lambda p: p**2/2 - gamma * np.log(p) + B*x - C, P_left))

print(np.allclose(res, p.compute_vertex_values(), rtol=1e-3))

plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
plot(p)

plt.plot(xs, res, ':r')
plt.legend(['numerical', 'exact'])

plt.show()
