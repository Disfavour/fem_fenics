"""
Solve for equation in 1D domain:

-a * laplacian(u) + b * (d(u)/dx) + c * u = f(x) + u(x)

with Essential BC on the sides
and Dirichlet BC u(1) = 0 on the right side
"""

from __future__ import print_function
from fenics import BoxMesh, Point, DirichletBC, FacetNormal
from fenics import FunctionSpace, VectorFunctionSpace
from fenics import TrialFunction, TestFunction
from fenics import Expression, Function, Constant, Measure
from fenics import dot, inner, Identity, sym
from fenics import grad, nabla_grad, tr
from ufl import nabla_div
from fenics import sqrt
from fenics import project
from fenics import dx,Dx
from fenics import solve, interpolate
from fenics import plot
from fenics import File
from fenics import lhs, rhs, assemble
from dolfin import *
import math
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('tkagg')


# Material settings
a = 1
b = 0.5
c = 2
pi = math.pi
alpha = 10e-2
gamma = alpha/2
# Domain settings
resW = 10  # mesh resolution (x, z)-directions (number of points)
print('Domain size:', 1, '[m] x')

# Finite element order
eo = 2

# Create a mesh
mesh = UnitIntervalMesh(resW) # it is an interval

# Create functional space
F = FunctionSpace(mesh, 'P', eo)  # piecewise linear polynomials

# Create 'trial function' (the unknown function to be approximated) over the space F
u = TrialFunction(F)  # scalar field, such that: u(x) = Sum(u_i* w(x))
us = TrialFunction(F)
manage = TrialFunction(F)
d1 = u.geometric_dimension()  # space dimension

# Create scalar 'test function' (or 'shape function', or 'weighting function') over the space F
q = TestFunction(F)  # q -> u

# Define expressions used in variational forms to prevent repeated code generation
n = FacetNormal(mesh)  # per-face normal vector
Traction = Constant(0)  # traction

source = '- cos(pi * x[0]) * (a*pi + (c/pi)) + ' \
        'sin(pi * x[0]) * b + ' \
        'x[0] * (2*a + b) +' \
        '(pow(x[0],2)) * ((c/2) - b)-' \
        '(pow(x[0], 3)) * c/3 -' \
        'a - c * ((1/pi) + (1/6))  '
source_f = Expression(source, degree=eo, a = a, b = b, c = c, pi = pi, eps=1.0e-9)

# Compute boundary surfaces of the mesh (via C++ wrapper)
tol = 1E-14

# exact solution
phi_obs = '- (1/pi)*cos(pi * x[0]) + (pow(x[0], 2))/2 - (pow(x[0], 3))/3 - ((1/pi) + (1/6))'
exact = Expression(phi_obs, degree=eo, pi = pi, eps=1.0e-9)
u_ex = interpolate(exact, F)

# Compose weak form of PDE(direct):
manage_0 = Expression('0', degree=eo) # manage equals zero
manage_prev = project(manage_0, F, solver_type="cg", preconditioner_type="amg")
u = TrialFunction(F)
q = TestFunction(F)

eq = (a * inner(grad(u), grad(q)) * dx)
eq += (b * Dx(u,0) * q * dx)
eq += (c * dot(u, q) * dx)
eq -= (dot(source_f + manage_prev, q) * dx)
a_u = lhs(eq)
L_u = rhs(eq)

# Compose weak form of PDE(undirect):
u_s = TrialFunction(F)
q = TestFunction(F)
eq_s = (a * inner(grad(u_s), grad(q)) * dx)
eq_s -= (b * Dx(u_s,0) * q * dx)
eq_s += (c * u_s * q * dx)
eq_s -= ((u - u_ex)* q * dx)
a_u_s = lhs(eq_s)
L_u_s = rhs(eq_s)

# Prepare export to file
path = 'results/'
u_file = File(f'{path}u.pvd')
u_ex_file = File(f'{path}u_ex.pvd')
manage_file = File(f'{path}manage.pvd')
u = Function(F)
u_s = Function(F)
manage = Function(F)
sol_settings = {'linear_solver': 'mumps'}

iters = 1
arr1 = []
while (iters < 1000):
        # Compute
        solve(a_u == L_u, u, solver_parameters=sol_settings)
        solve(a_u_s == L_u_s, u_s, solver_parameters=sol_settings)
       
        manage_f = 'manage_prev + gamma * (alpha * manage_prev + u_s) '
        manage_expression = Expression(manage_f, degree=eo, manage_prev = manage_prev, gamma = gamma, alpha = alpha, u_s = u_s, eps=1.0e-9)
        manage = interpolate(manage_expression, F)
        manage_array = manage.compute_vertex_values()
        manage_prev_array = manage_prev.compute_vertex_values()
        
        L2 = np.sqrt(((manage_array - manage_prev_array)**2).sum())
        if (iters != 1): 
                L2_u = np.sqrt(((manage_prev_array)**2).sum())
                if (iters%10 == 0):
                        gamma = gamma/2
                        arr1.append(L2/L2_u)

                if ((L2/L2_u) < 10e-3) and (iters!=1):
                        print(iters)
                        break
        manage_prev = manage
        print(iters)
        iters+= 1
# plt.plot([1,2,3])
print (arr1)
# Analyse
u_ex_array = u_ex.compute_vertex_values()
u_array = u.compute_vertex_values()

L2 = np.sqrt(((u_ex_array - u_array)**2).sum())
L1 = (np.abs(u_ex_array - u_array)).sum()

L2_u = np.sqrt(((u_ex_array)**2).sum())
L1_u = (np.abs(u_ex_array)).sum()

print('L1-error:', L1)
print('L1-related:', L1/L1_u)
print('L2-error:', L2)
print('L2-related:', L2/L2_u)


# Save solution
u.rename('u ', 'label')
u_file << (u)
u_ex.rename('u_ex ', 'label')
u_ex_file << (u_ex)
manage.rename('manage ', 'label')
manage_file << (manage)
