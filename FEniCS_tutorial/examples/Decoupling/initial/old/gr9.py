from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

plt.figure(num=0, figsize=(10, 6), dpi=160, facecolor='w', edgecolor='k')

tt = np.load("n50/K_1_0.01_t.npy")
y1 = np.load("n50/K_1_0.01_E.npy")
y2 = np.load("n50/K_5_0.01_E.npy")

plt.plot(tt, y1, 'C0--', label="$\\tau =0.01$")
plt.plot(tt, y2, 'C0-')

tt = np.load("n50/K_1_0.005_t.npy")
y1 = np.load("n50/K_1_0.005_E.npy")
y2 = np.load("n50/K_5_0.005_E.npy")

plt.plot(tt, y1, 'C1--', label="$\\tau =0.005$")
plt.plot(tt, y2, 'C1-')

plt.ylabel("$E$")
plt.xlabel("$t$")
plt.legend(loc=0)
plt.grid()

plt.show()

