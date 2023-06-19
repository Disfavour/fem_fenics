import matplotlib.pyplot as plt
import numpy as np

jac1t1 = np.load("../data/jac_K1_d1_ms100_tau0.01.npy")
jac2t1 = np.load("../data/jac_K2_d1_ms100_tau0.01.npy")
jac5t1 = np.load("../data/jac_K5_d1_ms100_tau0.01.npy")

jac1t2 = np.load("../data/jac_K1_d1_ms100_tau0.005.npy")
jac2t2 = np.load("../data/jac_K2_d1_ms100_tau0.005.npy")
jac5t2 = np.load("../data/jac_K5_d1_ms100_tau0.005.npy")

jac1t3 = np.load("../data/jac_K1_d1_ms100_tau0.0025.npy")
jac2t3 = np.load("../data/jac_K2_d1_ms100_tau0.0025.npy")
jac5t3 = np.load("../data/jac_K5_d1_ms100_tau0.0025.npy")

plt.figure(figsize=(16, 9), tight_layout=True)

plt.plot(jac1t1[:, 0], jac1t1[:, 1], color='blue', linestyle="dashed")
plt.plot(jac5t1[:, 0], jac5t1[:, 1], color='blue', linestyle="solid", label=r"$\tau=0.01$")

plt.plot(jac1t2[:, 0], jac1t2[:, 1], color='red', linestyle="dashed")
plt.plot(jac5t2[:, 0], jac5t2[:, 1], color='red', linestyle="solid", label=r"$\tau=0.005$")

plt.plot(jac1t3[:, 0], jac1t3[:, 1], color='green', linestyle="dashed")
plt.plot(jac5t3[:, 0], jac5t3[:, 1], color='green', linestyle="solid", label=r"$\tau=0.0025$")

plt.legend(loc='lower left')
plt.grid()
plt.xlim(jac1t1[:, 0].min(), jac1t1[:, 0].max())
plt.xlabel(r"$t$")
plt.ylabel(r"$E$")

plt.savefig('plots/fig11.pdf', dpi=300, bbox_inches='tight')
plt.show()
