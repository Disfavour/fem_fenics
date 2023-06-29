import matplotlib.pyplot as plt
import numpy as np


def my_plot(nl, jac1, jac2, jac5, filename):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    fig.set_size_inches(6.4, 3.6)
    fig.tight_layout()

    ax1.plot(nl[:, 0], nl[:, 1], color='blue', linestyle="solid")

    ax2.plot(nl[:, 0], np.abs(jac1[:, 1] - nl[:, 1]), color='red', linestyle="dotted", label=r"$K=1$")
    ax2.plot(nl[:, 0], np.abs(jac2[:, 1] - nl[:, 1]), color='red', linestyle="dashed", label=r"$K=2$")
    ax2.plot(nl[:, 0], np.abs(jac5[:, 1] - nl[:, 1]), color='red', linestyle="solid", label=r"$K=5$")

    ax2.set_yscale('log')

    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$E$')

    ax2.set_ylabel(r'$\delta E = |E^k - E|$')

    ax1.grid()
    #ax2.grid()

    #ax1.legend(loc='center left')
    #ax2.legend(loc='center right')
    #ax2.legend()

    ax1.set_xlim(nl[:, 0].min(), nl[:, 0].max())
    #ax2.set_ylim(ymin=1e-18)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')

    for i in [".pdf", ".png"]:
        plt.savefig(filename + i, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    for i in range(1, 4):
        nl = np.load(f"../data/nonlinear_d{i}_ms100_tau0.005.npy")

        jac1 = np.load(f"../data/jac_K1_d{i}_ms100_tau0.005.npy")
        jac2 = np.load(f"../data/jac_K2_d{i}_ms100_tau0.005.npy")
        jac5 = np.load(f"../data/jac_K5_d{i}_ms100_tau0.005.npy")

        sei1 = np.load(f"../data/sei_K1_d{i}_ms100_tau0.005.npy")
        sei2 = np.load(f"../data/sei_K2_d{i}_ms100_tau0.005.npy")
        sei5 = np.load(f"../data/sei_K5_d{i}_ms100_tau0.005.npy")

        my_plot(nl, jac1, jac2, jac5, f'plots/nl_jac_P{i}')
        my_plot(nl, sei1, sei2, sei5, f'plots/nl_sei_P{i}')

    for i in range(2, 4):
        nl = np.load(f"../data/nonlinear_RT_d{i}_ms100_tau0.005.npy")

        dec1 = np.load(f"../data/dec_RT_K1_d{i}_ms100_tau0.005.npy")
        dec2 = np.load(f"../data/dec_RT_K2_d{i}_ms100_tau0.005.npy")
        dec5 = np.load(f"../data/dec_RT_K5_d{i}_ms100_tau0.005.npy")

        my_plot(nl, dec1, dec2, dec5, f'plots/nl_dec_RT_{i}')
