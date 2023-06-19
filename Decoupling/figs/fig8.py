import matplotlib.pyplot as plt
import numpy as np


jac1 = np.load("../data/jac_K1_d1_ms100_tau0.005.npy")
jac2 = np.load("../data/jac_K2_d1_ms100_tau0.005.npy")
jac5 = np.load("../data/jac_K5_d1_ms100_tau0.005.npy")

plt.figure(figsize=(16, 9), tight_layout=True)

plt.plot(jac1[:, 2], jac1[:, 3], color='blue', linestyle="dashed")
plt.plot(jac1[:, 2], jac1[:, 4], color='red', linestyle="dashed")
plt.plot(jac1[:, 2], jac1[:, 5], color='green', linestyle="dashed")

plt.plot(jac2[:, 2], jac2[:, 3], color='blue', linestyle="dotted", linewidth=3)
plt.plot(jac2[:, 2], jac2[:, 4], color='red', linestyle="dotted", linewidth=3)
plt.plot(jac2[:, 2], jac2[:, 5], color='green', linestyle="dotted", linewidth=3)

plt.plot(jac5[:, 2], jac5[:, 3], color='blue', linestyle="solid", label=r"$t=1$")
plt.plot(jac5[:, 2], jac5[:, 4], color='red', linestyle="solid", label=r"$t=2$")
plt.plot(jac5[:, 2], jac5[:, 5], color='green', linestyle="solid", label=r"$t=3$")

plt.legend()
plt.grid()
plt.xlim(jac1[:, 2].min(), jac1[:, 2].max())
plt.xlabel(r"$x_1$")
plt.ylabel(r"$\varrho$")

plt.savefig('plots/fig8.pdf', dpi=300, bbox_inches='tight')
plt.show()
