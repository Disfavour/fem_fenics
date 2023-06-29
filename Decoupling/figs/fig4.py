import matplotlib.pyplot as plt
import numpy as np


def fig4(nl_t1, nl_t2, nl_t3, filename):
    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

    plt.plot(nl_t1[:, 2], nl_t1[:, 3], color='blue', linestyle="dotted")
    plt.plot(nl_t1[:, 2], nl_t1[:, 4], color='red', linestyle="dotted")
    plt.plot(nl_t1[:, 2], nl_t1[:, 5], color='green', linestyle="dotted")

    plt.plot(nl_t2[:, 2], nl_t2[:, 3], color='blue', linestyle="dashed")
    plt.plot(nl_t2[:, 2], nl_t2[:, 4], color='red', linestyle="dashed")
    plt.plot(nl_t2[:, 2], nl_t2[:, 5], color='green', linestyle="dashed")

    plt.plot(nl_t3[:, 2], nl_t3[:, 3], color='blue', linestyle="solid", label=r"$t=1$")
    plt.plot(nl_t3[:, 2], nl_t3[:, 4], color='red', linestyle="solid", label=r"$t=2$")
    plt.plot(nl_t3[:, 2], nl_t3[:, 5], color='green', linestyle="solid", label=r"$t=3$")

    plt.legend()
    plt.grid()
    plt.xlim(nl_t1[:, 2].min(), nl_t1[:, 2].max())
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$\varrho$")

    for i in [".pdf", ".png"]:
        plt.savefig(filename + i, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    for i in range(1, 4):
        nl_t1 = np.load(f"../data/nonlinear_d{i}_ms100_tau0.01.npy")
        nl_t2 = np.load(f"../data/nonlinear_d{i}_ms100_tau0.005.npy")
        nl_t3 = np.load(f"../data/nonlinear_d{i}_ms100_tau0.0025.npy")

        fig4(nl_t1, nl_t2, nl_t3, f'plots/fig4_P{i}')

    for i in range(2, 4):
        nl_t1 = np.load(f"../data/nonlinear_RT_d{i}_ms100_tau0.01.npy")
        nl_t2 = np.load(f"../data/nonlinear_RT_d{i}_ms100_tau0.005.npy")
        nl_t3 = np.load(f"../data/nonlinear_RT_d{i}_ms100_tau0.0025.npy")

        fig4(nl_t1, nl_t2, nl_t3, f'plots/fig4_RT_{i}')
