from fenics import *
import numpy as np
import matplotlib.pyplot as plt


Theta = pi / 2
a, b = 1, 5.0
nr = 20
nt = 40
mesh = RectangleMesh(Point(a, 0), Point(b, 1), nr, nt, "right")


x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
s = 1.5

def denser(x, y):
    return [a + (b - a) * ((x - a) / (b - a)) ** s, y]

x_bar, y_bar = denser(x, y)
xy_bar_coor = np.array([x_bar, y_bar]).transpose()
mesh.coordinates()[:] = xy_bar_coor

plot(mesh)
plt.show()
print(mesh.num_cells())
print(refine(mesh).num_cells())
print(refine(refine(mesh)).num_cells())

def cylinder(r, s):
    return [r*np.cos(Theta*s), r*np.sin(Theta*s)]

x_hat, y_hat = cylinder(x_bar, y_bar)
xy_hat_coor = np.array([x_hat, y_hat]).transpose()
mesh.coordinates()[:] = xy_hat_coor

plot(mesh)
#plt.savefig("3.png")
plt.show()
print(mesh.num_cells())

mesh = refine(mesh)
print(mesh.num_cells())

plot(mesh)
#plt.savefig("2.png")
plt.show()
