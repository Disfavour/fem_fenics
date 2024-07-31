import pipe, pipe_linearized
import matplotlib.pyplot as plt
from os.path import join
import numpy as np

taus = [i * 3600 / 2 for i in (1, 2, 4)]
print(taus)

lines = ['-', '--', ':', '-.']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

ke_ti = np.genfromtxt(join('data', 'pipes', 'periodic_Ke-Ti.csv'), delimiter=',')
Osiadacz = np.genfromtxt(join('data', 'pipes', 'periodic_Osiadacz.csv'), delimiter=',')

plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

tau = taus[0]
t, P_in, P_out, m_in, m_out = pipe.calculate_pipe(mesh_size=100, tau=tau)
plt.plot(t, P_out, colors[0], label='implicit')
t, P_in, P_out, m_in, m_out = pipe_linearized.calculate_pipe(mesh_size=100, tau=tau)
plt.plot(t, P_out, colors[1], label='linearized')

for tau, l in zip(taus[1:], lines[1:]):
    t, P_in, P_out, m_in, m_out = pipe.calculate_pipe(mesh_size=100, tau=tau)
    plt.plot(t, P_out, colors[0]+l)
    t, P_in, P_out, m_in, m_out = pipe_linearized.calculate_pipe(mesh_size=100, tau=tau)
    plt.plot(t, P_out, colors[1]+l)

plt.plot(ke_ti[:,0], ke_ti[:,1], 'o'+colors[2], ms=3, label="Ke & Ti")
plt.plot(Osiadacz[:,0], Osiadacz[:,1], 'o'+colors[3], ms=3, label='Osiadacz')

plt.xticks(range(0, 49*3600, 8*3600))
plt.xlabel(r'$t$')
plt.ylabel(r'$p_{out}$')
plt.xlim(t[0], t[-1])
plt.legend()
plt.grid()
plt.savefig(join('images', 'pipes', 'periodic.pdf'), transparent=True)
