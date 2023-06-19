from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

plt.figure(num=0, figsize=(10, 6), dpi=160, facecolor='w', edgecolor='k')

x = np.load("ref_x.npy")
y1 = np.load("tau_0.01_1.npy")
y2 = np.load("tau_0.01_2.npy")
y3 = np.load("tau_0.01_3.npy")

plt.plot(x, y1, 'C0:')
plt.plot(x, y2, 'C1:')
plt.plot(x, y3, 'C2:')

y1 = np.load("tau_0.005_1.npy")
y2 = np.load("tau_0.005_2.npy")
y3 = np.load("tau_0.005_3.npy")

plt.plot(x, y1, 'C0--')
plt.plot(x, y2, 'C1--')
plt.plot(x, y3, 'C2--')

y1 = np.load("tau_0.0025_1.npy")
y2 = np.load("tau_0.0025_2.npy")
y3 = np.load("tau_0.0025_3.npy")

plt.plot(x, y1, 'C0', label="$t = 1$")
plt.plot(x, y2, 'C1', label="$t = 2$")
plt.plot(x, y3, 'C2', label="$t = 3$")

plt.ylabel("$\\varrho$")
plt.xlabel("$x_1$")  
plt.legend(loc=0)
plt.ylim(0.86,1.24)
plt.grid()

plt.show()

