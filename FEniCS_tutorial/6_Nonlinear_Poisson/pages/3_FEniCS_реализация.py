import streamlit as st


r"""
# FEniCS реализация

```python
from fenics import *
import matplotlib.pyplot as plt


u_e = Expression("1 + x[0] + 2 * x[1]", degree=2)

mesh = UnitSquareMesh(8, 8, "right")
V = FunctionSpace(mesh, 'P', 1)

bc = DirichletBC(V, u_e, "on_boundary")

u = Function(V)
v = TestFunction(V)
f = Expression("-10 - 10*x[0] - 20*x[1]", degree=2)


def q(u):
    return 1 + u**2


F = (q(u) * dot(grad(u), grad(v)) - f * v) * dx

solve(F == 0, u, bc)

error_L2 = errornorm(u_e, u)

plot(u)
plt.show()
```
"""
