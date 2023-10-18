import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.ticker as mtick


base_dir = os.path.dirname(os.path.dirname(__file__))
data = os.path.join(base_dir, 'data')
plots = os.path.join(base_dir, 'plots')

fig, axs = plt.subplots(2, 2)

#fig.set_size_inches(6.4, 3.6)
fig.set_size_inches(10, 6)
#fig.tight_layout()
#fig.set_dpi(600)

lws = [1.5, 1.2, 0.9]


s1_t1 = np.load(os.path.join(data, f'w_s1_tau0.01_ms200.npy'))
s1_t2 = np.load(os.path.join(data, f'w_s1_tau0.005_ms200.npy'))
s1_t3 = np.load(os.path.join(data, f'w_s1_tau0.0025_ms200.npy'))
s2_t1 = np.load(os.path.join(data, f'w_s0.5_tau0.01_ms200.npy'))
s2_t2 = np.load(os.path.join(data, f'w_s0.5_tau0.005_ms200.npy'))
s2_t3 = np.load(os.path.join(data, f'w_s0.5_tau0.0025_ms200.npy'))

axs[0, 0].plot(s1_t1[10], s1_t1[11], ':b', lw=lws[0])
axs[0, 0].plot(s1_t1[10], s1_t1[12], ':r', lw=lws[0])
axs[0, 0].plot(s1_t1[10], s1_t1[13], ':g', lw=lws[0])

axs[0, 0].plot(s1_t2[10], s1_t2[11], '--b', lw=lws[1])
axs[0, 0].plot(s1_t2[10], s1_t2[12], '--r', lw=lws[1])
axs[0, 0].plot(s1_t2[10], s1_t2[13], '--g', lw=lws[1])

axs[0, 0].plot(s1_t3[10], s1_t3[11], '-b', lw=lws[2], label=r"$t=1$")
axs[0, 0].plot(s1_t3[10], s1_t3[12], '-r', lw=lws[2], label=r"$t=2$")
axs[0, 0].plot(s1_t3[10], s1_t3[13], '-g', lw=lws[2], label=r"$t=3$")

axs[0, 0].set_xlim(s1_t1[10][0], s1_t1[10][-1])
axs[0, 0].legend()
axs[0, 0].grid()
axs[0, 0].set_xlabel(r'$x_1$')
axs[0, 0].set_ylabel(r'$\varrho$')

axs[0, 1].plot(s2_t1[10], s2_t1[11], ':b', lw=lws[0])
axs[0, 1].plot(s2_t1[10], s2_t1[12], ':r', lw=lws[0])
axs[0, 1].plot(s2_t1[10], s2_t1[13], ':g', lw=lws[0])

axs[0, 1].plot(s2_t2[10], s2_t2[11], '--b', lw=lws[1])
axs[0, 1].plot(s2_t2[10], s2_t2[12], '--r', lw=lws[1])
axs[0, 1].plot(s2_t2[10], s2_t2[13], '--g', lw=lws[1])

axs[0, 1].plot(s2_t3[10], s2_t3[11], '-b', lw=lws[2], label=r"$t=1$")
axs[0, 1].plot(s2_t3[10], s2_t3[12], '-r', lw=lws[2], label=r"$t=2$")
axs[0, 1].plot(s2_t3[10], s2_t3[13], '-g', lw=lws[2], label=r"$t=3$")

axs[0, 1].set_xlim(s1_t1[10][0], s1_t1[10][-1])
axs[0, 1].legend()
axs[0, 1].grid()
axs[0, 1].set_xlabel(r'$x_1$')
axs[0, 1].set_ylabel(r'$\varrho$')

axs[1, 0].plot(s1_t1[0], s1_t1[2], ':b', lw=1.5, label=r'$\tau=0.01$')
axs[1, 0].plot(s1_t2[0], s1_t2[2], '--r', lw=1.2, label=r"$\tau=0.005$")
axs[1, 0].plot(s1_t3[0], s1_t3[2], '-g', lw=0.9, label=r"$\tau=0.0025$")

axs[1, 0].set_xlim(s1_t1[0][0], s1_t1[0][-1])
axs[1, 0].legend()
axs[1, 0].grid()
axs[1, 0].set_xlabel(r'$t$')
axs[1, 0].set_ylabel(r'$E$')

axs[1, 1].plot(s2_t1[0], s2_t1[2], ':b', lw=1.5, label=r'$\tau=0.01$')
axs[1, 1].plot(s2_t2[0], s2_t2[2], '--r', lw=1.2, label=r"$\tau=0.005$")
axs[1, 1].plot(s2_t3[0], s2_t3[2], '-g', lw=0.9, label=r"$\tau=0.0025$")

axs[1, 1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))

axs[1, 1].set_xlim(s2_t1[0][0], s2_t1[0][-1])
axs[1, 1].legend()
axs[1, 1].grid()
axs[1, 1].set_xlabel(r'$t$')
axs[1, 1].set_ylabel(r'$E$')
fig.tight_layout()
fig.savefig(os.path.join(plots, 'tihonov.png'), dpi=300)

fig.show()

# if arr_t1[2].max() - arr_t1[2].min() < 0.05:
#     plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))


