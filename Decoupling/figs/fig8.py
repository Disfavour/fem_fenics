import matplotlib.pyplot as plt
import numpy as np


def fig8(jac1, jac2, jac5, filename):
    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

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

    for i in [".pdf", ".png"]:
        plt.savefig(filename + i, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    for i in range(1, 4):
        jac1 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.005.npy")
        jac2 = np.load(f"../data/jac_K2_d{i}_ms100_tau0.005.npy")
        jac5 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.005.npy")

        fig8(jac1, jac2, jac5, f'plots/fig8_P{i}')

    for i in range(1, 4):
        dec_RT_K1 = np.load(f"../data/dec_RT_K1_d{i}_ms100_tau0.005.npy")
        dec_RT_K2 = np.load(f"../data/dec_RT_K2_d{i}_ms100_tau0.005.npy")
        dec_RT_K5 = np.load(f"../data/dec_RT_K5_d{i}_ms100_tau0.005.npy")

        fig8(dec_RT_K1, dec_RT_K2, dec_RT_K5, f'plots/fig8_RT_{i}')
