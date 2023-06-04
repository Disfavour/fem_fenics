from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

plt.figure(num=0, figsize=(10, 6), dpi=160, facecolor='w', edgecolor='k')

x = np.load("ref_x.npy")
y1 = np.load("K_1_1.npy")
y2 = np.load("K_1_2.npy")
y3 = np.load("K_1_3.npy")

plt.plot(x, y1, 'C0--')
plt.plot(x, y2, 'C1--')
plt.plot(x, y3, 'C2--')

y1 = np.load("K_2_1.npy")
y2 = np.load("K_2_2.npy")
y3 = np.load("K_2_3.npy")

plt.plot(x, y1, 'C0.')
plt.plot(x, y2, 'C1.')
plt.plot(x, y3, 'C2.')

y1 = np.load("K_5_1.npy")
y2 = np.load("K_5_2.npy")
y3 = np.load("K_5_3.npy")

plt.plot(x, y1, 'C0', label="$t = 1$")
plt.plot(x, y2, 'C1', label="$t = 2$")
plt.plot(x, y3, 'C2', label="$t = 3$")

plt.ylabel("$\\varrho$")
plt.xlabel("$x_1$")  
plt.legend(loc=0)
plt.ylim(0.86,1.24)
plt.grid()

plt.show()

