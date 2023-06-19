import matplotlib.pyplot as plt
import numpy as np


nl = np.load("../data/nonlinear_d1_ms100_tau0.005.npy")
sei1 = np.load("../data/sei_K1_d1_ms100_tau0.005.npy")
sei2 = np.load("../data/sei_K2_d1_ms100_tau0.005.npy")
sei5 = np.load("../data/sei_K5_d1_ms100_tau0.005.npy")


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

fig.set_size_inches(16, 9)
fig.tight_layout()

ax1.plot(nl[:, 0], nl[:, 1], color='blue', linestyle="solid", label=r"nonlinear")

ax2.plot(nl[:, 0], sei1[:, 1] - nl[:, 1], color='red', linestyle="dashed", label=r"$K=1$")
ax2.plot(nl[:, 0], sei2[:, 1] - nl[:, 1], color='red', linestyle="dotted", label=r"$K=2$", linewidth=3)
ax2.plot(nl[:, 0], sei5[:, 1] - nl[:, 1], color='red', linestyle="solid", label=r"$K=5$")

ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$E$')

ax2.set_ylabel(r'$\partial E = E^k - E$')

ax1.grid()
#ax2.grid()

ax1.legend(loc='center left')
ax2.legend(loc='center right')

ax1.set_xlim(nl[:, 0].min(), nl[:, 0].max())

plt.savefig('plots/nl_sei.pdf', dpi=300, bbox_inches='tight')
plt.show()
