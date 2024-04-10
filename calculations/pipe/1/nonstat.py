from fenics import *
import numpy as np
from scipy.constants import pi, R
import matplotlib.pyplot as plt
from  scipy.optimize import newton


#set_log_level(LogLevel.WARNING)

P_in, P_out, m_in, m_out = [], [], [], []

t_max = 1
mesh_size = 100
tau = 0.1

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
m_right = rho * 70


ts = np.arange(0, t_max+tau/2, tau)
mesh = IntervalMesh(mesh_size, 0, L)

gamma = Rs * T * m_right**2 / A**2
B = f * m_right**2 * Rs * T / (2 * D * A**2)
C = P_left**2 / 2 - gamma * np.log(P_left)

xs = mesh.coordinates().flatten()
exact = []
for x in xs:
    exact.append(newton(lambda p: p**2/2 - gamma * np.log(p) + B*x - C, P_left))
exact = np.array(exact)

P = FiniteElement('P', mesh.ufl_cell(), 1)
M = FiniteElement('P', mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, MixedElement([P, M]))

bc = [DirichletBC(W.sub(0), Constant(P_left), 'on_boundary && near(x[0], 0)'),
      DirichletBC(W.sub(1), Constant(m_right), CompiledSubDomain('on_boundary && near(x[0], L)', L=L))]

w, wn = Function(W), Function(W)
pt, mt = TestFunctions(W)
p, m = split(w)
pn, mn = split(wn)

F = ((p-pn)/(tau*Rs*T) + m.dx(0)/A) * pt*dx \
    + ((m-mn)/(tau*A) + Rs*T/A**2*(m**2/p).dx(0) + p.dx(0) + Rs*T*f*m*abs(m)/(2*D*A**2*p)) * mt*dx


t = 0
w.assign(project(Expression(('P_left-30*x[0]', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W))
#w.assign(project(Expression(('P_left', 'm_right'), P_left=P_left, m_right=m_right, degree=1), W))

#w.vector()[dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()]] = exact
wn.assign(w)

# plt.figure()
# plot(w.sub(0))
# plt.figure()
# plot(w.sub(1))
# plt.show()
# exit()




for t in ts[1:]:
    solve(F == 0, w, bc, solver_parameters={"newton_solver": {
            'absolute_tolerance': 1e-5,
            'relative_tolerance': 1e-8,
            'maximum_iterations': 500,
            'relaxation_parameter': 1.0,
        }})
    wn.assign(w)

print(np.allclose(exact, w.sub(0).compute_vertex_values(), rtol=1e-3))

plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
plot(w.sub(0))

plt.plot(xs, exact, ':r')
plt.legend(['numerical', 'exact'])

plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)
plot(w.sub(1))

plt.show()
