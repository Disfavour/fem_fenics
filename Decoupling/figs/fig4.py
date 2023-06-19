import matplotlib.pyplot as plt
import numpy as np


nl1 = np.load("../data/nonlinear_d1_ms100_tau0.01.npy")
nl2 = np.load("../data/nonlinear_d1_ms100_tau0.005.npy")
nl3 = np.load("../data/nonlinear_d1_ms100_tau0.0025.npy")

plt.figure(figsize=(16, 9), tight_layout=True)

plt.plot(nl1[:, 2], nl1[:, 3], color='blue', linestyle="dotted")
plt.plot(nl1[:, 2], nl1[:, 4], color='red', linestyle="dotted")
plt.plot(nl1[:, 2], nl1[:, 5], color='green', linestyle="dotted")

plt.plot(nl2[:, 2], nl2[:, 3], color='blue', linestyle="dashed")
plt.plot(nl2[:, 2], nl2[:, 4], color='red', linestyle="dashed")
plt.plot(nl2[:, 2], nl2[:, 5], color='green', linestyle="dashed")

plt.plot(nl3[:, 2], nl3[:, 3], color='blue', linestyle="solid", label=r"$t=1$")
plt.plot(nl3[:, 2], nl3[:, 4], color='red', linestyle="solid", label=r"$t=2$")
plt.plot(nl3[:, 2], nl3[:, 5], color='green', linestyle="solid", label=r"$t=3$")

plt.legend()
plt.grid()
plt.xlim(nl3[:, 2].min(), nl3[:, 2].max())
plt.xlabel(r"$x_1$")
plt.ylabel(r"$\varrho$")

plt.savefig('plots/fig4.pdf', dpi=300, bbox_inches='tight')
plt.show()
