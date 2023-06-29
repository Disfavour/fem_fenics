import matplotlib.pyplot as plt
import numpy as np


def fig11(jac1_t1, jac5_t1, jac1_t2, jac5_t2, jac1_t3, jac5_t3, filename):
    plt.figure(figsize=(6.4, 3.6), tight_layout=True)

    plt.plot(jac1_t1[:, 0], jac1_t1[:, 1], color='blue', linestyle="dashed")
    plt.plot(jac5_t1[:, 0], jac5_t1[:, 1], color='blue', linestyle="solid", label=r"$\tau=0.01$")

    plt.plot(jac1_t2[:, 0], jac1_t2[:, 1], color='red', linestyle="dashed")
    plt.plot(jac5_t2[:, 0], jac5_t2[:, 1], color='red', linestyle="solid", label=r"$\tau=0.005$")

    plt.plot(jac1_t3[:, 0], jac1_t3[:, 1], color='green', linestyle="dashed")
    plt.plot(jac5_t3[:, 0], jac5_t3[:, 1], color='green', linestyle="solid", label=r"$\tau=0.0025$")

    plt.legend()
    plt.grid()
    plt.xlim(jac1_t1[:, 0].min(), jac1_t1[:, 0].max())
    plt.xlabel(r"$t$")
    plt.ylabel(r"$E$")

    for i in [".pdf", ".png"]:
        plt.savefig(filename + i, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    for i in range(1, 4):
        jac1_t1 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.01.npy")
        jac5_t1 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.01.npy")

        jac1_t2 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.005.npy")
        jac5_t2 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.005.npy")

        jac1_t3 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.0025.npy")
        jac5_t3 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.0025.npy")

        fig11(jac1_t1, jac5_t1, jac1_t2, jac5_t2, jac1_t3, jac5_t3, f'plots/fig11_P{i}')

    for i in range(1, 4):
        dec_RT_K1_t1 = np.load(f"../data/dec_RT_K1_d{i}_ms100_tau0.01.npy")
        dec_RT_K5_t1 = np.load(f"../data/dec_RT_K5_d{i}_ms100_tau0.01.npy")

        dec_RT_K1_t2 = np.load(f"../data/dec_RT_K1_d{i}_ms100_tau0.005.npy")
        dec_RT_K5_t2 = np.load(f"../data/dec_RT_K5_d{i}_ms100_tau0.005.npy")

        dec_RT_K1_t3 = np.load(f"../data/dec_RT_K1_d{i}_ms100_tau0.0025.npy")
        dec_RT_K5_t3 = np.load(f"../data/dec_RT_K5_d{i}_ms100_tau0.0025.npy")

        fig11(dec_RT_K1_t1, dec_RT_K5_t1, dec_RT_K1_t2, dec_RT_K5_t2, dec_RT_K1_t3, dec_RT_K5_t3, f'plots/fig11_RT_{i}')
