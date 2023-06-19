import matplotlib.pyplot as plt
import numpy as np


nl = np.load("../data/nonlinear_d1_ms100_tau0.005.npy")
jac1 = np.load("../data/jac_K1_d1_ms100_tau0.005.npy")
jac2 = np.load("../data/jac_K2_d1_ms100_tau0.005.npy")
jac5 = np.load("../data/jac_K5_d1_ms100_tau0.005.npy")

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

fig.set_size_inches(16, 9)
fig.tight_layout()

ax1.plot(nl[:, 0], nl[:, 1], color='blue', linestyle="solid", label=r"nonlinear")

ax2.plot(nl[:, 0], jac1[:, 1] - nl[:, 1], color='red', linestyle="dashed", label=r"$K=1$")
ax2.plot(nl[:, 0], jac2[:, 1] - nl[:, 1], color='red', linestyle="dotted", label=r"$K=2$", linewidth=3)
ax2.plot(nl[:, 0], jac5[:, 1] - nl[:, 1], color='red', linestyle="solid", label=r"$K=5$")

ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$E$')

ax2.set_ylabel(r'$\partial E = E^k - E$')

ax1.grid()
#ax2.grid()

ax1.legend(loc='center left')
ax2.legend(loc='center right')

ax1.set_xlim(nl[:, 0].min(), nl[:, 0].max())

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
#ax1.legend(lines1 + lines2, labels1 + labels2)

plt.savefig('plots/nl_jac.pdf', dpi=300, bbox_inches='tight')
plt.show()
