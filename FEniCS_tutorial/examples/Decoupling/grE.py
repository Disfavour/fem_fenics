from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

plt.figure(num=0, figsize=(10, 6), dpi=160, facecolor='w', edgecolor='k')

tt = np.load("K_1_0.01_t.npy")
y1 = np.load("K_1_0.01_E.npy")
y2 = np.load("K_5_0.01_E.npy")

plt.plot(tt, y1, 'C0--')
plt.plot(tt, y2, 'C0-', label="$\\tau =0.01$")

tt = np.load("K_1_0.005_t.npy")
y1 = np.load("K_1_0.005_E.npy")
y2 = np.load("K_5_0.005_E.npy")

plt.plot(tt, y1, 'C1--')
plt.plot(tt, y2, 'C1-', label="$\\tau =0.005$")

tt = np.load("K_1_0.0025_t.npy")
y1 = np.load("K_1_0.0025_E.npy")
y2 = np.load("K_5_0.0025_E.npy")

plt.plot(tt, y1, 'C2--')
plt.plot(tt, y2, 'C2-', label="$\\tau =0.0025$")

plt.ylabel("$E$")
plt.xlabel("$t$")
plt.legend(loc=3)
plt.grid()
plt.ylim(251.1525,251.2875)

plt.show()

