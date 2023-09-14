import streamlit as st
import os.path
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt


files = os.path.join(os.path.dirname(os.path.dirname(__file__)), "files")

vis = open(os.path.join(files, "vis.mp4"), 'rb').read()
im = Image.open(os.path.join(files, "2.png"))


def fig4(nl_t1, nl_t2, nl_t3):
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

    return plt.gcf()


def fig8(jac1, jac2, jac5):
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

    return plt.gcf()


def fig11(jac1_t1, jac5_t1, jac1_t2, jac5_t2, jac1_t3, jac5_t3):
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

    return plt.gcf()


def my_plot(nl, jac1, jac2, jac5):
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

    ax1.set_xlim(nl[:, 0].min(), nl[:, 0].max())

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')

    return plt.gcf()


menu = st.sidebar.radio(
    '***',
    (
        'Динамика плотности',
        'Динамика плотностей',
        'Плотность на различные моменты времени (Неявная схема)',
        'Плотность на различные моменты времени (Итерационная cхема)',
        'Полная механическая энергия на различные моменты времени',
        'Отклонение полной механической энергии',
    )
)

if menu == 'Динамика плотности':
    r"""
    ## Динамика плотности
    """
    st.video(vis)

if menu == 'Динамика плотностей':
    r"""
    ## Динамика плотности в центре, максимального и минимального значений плотности
    """
    st.image(im)

if menu == 'Плотность на различные моменты времени (Неявная схема)':
    r"""
    ## Плотность на различные моменты времени (Неявная схема)
    """
    # nl_t1 = np.load(os.path.join(files, 'nonlinear_d1_ms100_tau0.01.npy'))
    # nl_t2 = np.load(os.path.join(files, 'nonlinear_d1_ms100_tau0.005.npy'))
    # nl_t3 = np.load(os.path.join(files, 'nonlinear_d1_ms100_tau0.0025.npy'))
    # st.pyplot(fig4(nl_t1, nl_t2, nl_t3))
    r'Линия из точек — $\tau = 0.01$ , штриховая линия — $\tau = 0.005$, сплошная — $\tau = 0.0025$'

    for i in range(1, 4):
        f'##### Лагранжевы элементы {i}-й степени'
        nl_t1 = np.load(os.path.join(files, f'nonlinear_d{i}_ms100_tau0.01.npy'))
        nl_t2 = np.load(os.path.join(files, f'nonlinear_d{i}_ms100_tau0.005.npy'))
        nl_t3 = np.load(os.path.join(files, f'nonlinear_d{i}_ms100_tau0.0025.npy'))
        st.pyplot(fig4(nl_t1, nl_t2, nl_t3))

if menu == 'Плотность на различные моменты времени (Итерационная cхема)':
    r"""
    ## Плотность на различные моменты времени (Итерационная cхема)
    """
    # jac1 = np.load(os.path.join(files, 'jac_K1_d1_ms100_tau0.005.npy'))
    # jac2 = np.load(os.path.join(files, 'jac_K2_d1_ms100_tau0.005.npy'))
    # jac5 = np.load(os.path.join(files, 'jac_K5_d1_ms100_tau0.005.npy'))
    # st.pyplot(fig8(jac1, jac2, jac5))
    r'Штриховая линия — $K = 1$, линия из точек — $K = 2$, сплошная — $K = 5$'

    for i in range(1, 4):
        f'##### Лагранжевы элементы {i}-й степени'
        jac1 = np.load(os.path.join(files, f'jac_K1_d{i}_ms100_tau0.005.npy'))
        jac2 = np.load(os.path.join(files, f'jac_K2_d{i}_ms100_tau0.005.npy'))
        jac5 = np.load(os.path.join(files, f'jac_K5_d{i}_ms100_tau0.005.npy'))
        st.pyplot(fig8(jac1, jac2, jac5))

if menu == 'Полная механическая энергия на различные моменты времени':
    r"""
    ## Полная механическая энергия на различные моменты времени
    """
    # jac1_t1 = np.load(os.path.join(files, 'jac_K1_d1_ms100_tau0.01.npy'))
    # jac5_t1 = np.load(os.path.join(files, 'jac_K5_d1_ms100_tau0.01.npy'))
    # jac1_t2 = np.load(os.path.join(files, 'jac_K1_d1_ms100_tau0.005.npy'))
    # jac5_t2 = np.load(os.path.join(files, 'jac_K5_d1_ms100_tau0.005.npy'))
    # jac1_t3 = np.load(os.path.join(files, 'jac_K1_d1_ms100_tau0.0025.npy'))
    # jac5_t3 = np.load(os.path.join(files, 'jac_K5_d1_ms100_tau0.0025.npy'))
    # st.pyplot(fig11(jac1_t1, jac5_t1, jac1_t2, jac5_t2, jac1_t3, jac5_t3))
    r'Штриховая линия — $K = 1$, сплошная — $K = 5$'

    for i in range(1, 4):
        f'##### Лагранжевы элементы {i}-й степени'
        jac1_t1 = np.load(os.path.join(files, f'jac_K1_d{i}_ms100_tau0.01.npy'))
        jac5_t1 = np.load(os.path.join(files, f'jac_K5_d{i}_ms100_tau0.01.npy'))
        jac1_t2 = np.load(os.path.join(files, f'jac_K1_d{i}_ms100_tau0.005.npy'))
        jac5_t2 = np.load(os.path.join(files, f'jac_K5_d{i}_ms100_tau0.005.npy'))
        jac1_t3 = np.load(os.path.join(files, f'jac_K1_d{i}_ms100_tau0.0025.npy'))
        jac5_t3 = np.load(os.path.join(files, f'jac_K5_d{i}_ms100_tau0.0025.npy'))
        st.pyplot(fig11(jac1_t1, jac5_t1, jac1_t2, jac5_t2, jac1_t3, jac5_t3))

if menu == 'Отклонение полной механической энергии':
    r"""
    ## Отклонение полной механической энергии
    """
    nl = np.load(os.path.join(files, 'nonlinear_d1_ms100_tau0.005.npy'))
    jac1 = np.load(os.path.join(files, 'jac_K1_d1_ms100_tau0.005.npy'))
    jac2 = np.load(os.path.join(files, 'jac_K2_d1_ms100_tau0.005.npy'))
    jac5 = np.load(os.path.join(files, 'jac_K5_d1_ms100_tau0.005.npy'))
    sei1 = np.load(os.path.join(files, 'sei_K1_d1_ms100_tau0.005.npy'))
    sei2 = np.load(os.path.join(files, 'sei_K2_d1_ms100_tau0.005.npy'))
    sei5 = np.load(os.path.join(files, 'sei_K5_d1_ms100_tau0.005.npy'))
    '### Расщепление типа Якоби'
    st.pyplot(my_plot(nl, jac1, jac2, jac5))
    '### Расщепление типа Зейделя'
    st.pyplot(my_plot(nl, sei1, sei2, sei5))
