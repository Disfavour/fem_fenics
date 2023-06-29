import matplotlib.pyplot as plt
import numpy as np


def dif(nl, jac1, jac2, jac5, sei1, sei2, sei5, filename):
    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

    plt.plot(jac1[:, 0], np.abs(jac1[:, 1] - nl[:, 1]), color='blue', linestyle="dashdot", label=r'jac $K=1$')
    plt.plot(jac2[:, 0], np.abs(jac2[:, 1] - nl[:, 1]), color='blue', linestyle="dashed", label=r'jac $K=2$')
    plt.plot(jac5[:, 0], np.abs(jac5[:, 1] - nl[:, 1]), color='blue', linestyle="solid", label=r'jac $K=5$')

    plt.plot(sei1[:, 0], np.abs(sei1[:, 1] - nl[:, 1]), color='red', linestyle="dotted", label=r'sei $K=1$', linewidth=3)
    plt.plot(sei2[:, 0], np.abs(sei2[:, 1] - nl[:, 1]), color='red', linestyle="dashed", label=r'sei $K=2$')
    plt.plot(sei5[:, 0], np.abs(sei5[:, 1] - nl[:, 1]), color='red', linestyle="solid", label=r'sei $K=5$')

    plt.yscale('log')

    plt.legend()
    plt.grid()
    plt.xlim(nl[:, 0].min(), nl[:, 0].max())
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\delta E = |E^k - E|$")

    plt.savefig(filename, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    for i in range(1, 4):
        nl = np.load(f"../data/nonlinear_d{i}_ms100_tau0.005.npy")

        jac1 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.005.npy")
        jac2 = np.load(f"../data/jac_K2_d{i}_ms100_tau0.005.npy")
        jac5 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.005.npy")

        sei1 = np.load(f"../data/sei_K1_d{i}_ms100_tau0.005.npy")
        sei2 = np.load(f"../data/sei_K2_d{i}_ms100_tau0.005.npy")
        sei5 = np.load(f"../data/sei_K5_d{i}_ms100_tau0.005.npy")

        dif(nl, jac1, jac2, jac5, sei1, sei2, sei5, f'plots/dif_P{i}.pdf')
