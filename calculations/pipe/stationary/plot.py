import nonstationary, nonstationary_linearized
import matplotlib.pyplot as plt
from os.path import join

taus = [i * 3600 / 16 for i in (1, 4, 16)]
tmax=taus[-1] * 2
print(taus, tmax)

lines = ['-', '--', ':', '-.']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

plt.figure(figsize=(6.4, 3.6), dpi=300, tight_layout=True)

x, p_numerical, p_analytic = nonstationary.calculate(mesh_size=100, tau=taus[0], t_max=tmax)
plt.plot(x, p_numerical, colors[0], label="implicit")
x, p_numerical, p_analytic = nonstationary_linearized.calculate(mesh_size=100, tau=taus[0], t_max=tmax)
plt.plot(x, p_numerical, colors[1], label="linearized")
plt.plot(x, p_analytic, colors[2], label='exact')

for tau, l in zip(taus[1:], lines[1:]):
    x, p_numerical, p_analytic = nonstationary.calculate(mesh_size=100, tau=tau, t_max=tmax)
    plt.plot(x, p_numerical, colors[0]+l)
    x, p_numerical, p_analytic = nonstationary_linearized.calculate(mesh_size=100, tau=tau, t_max=tmax)
    plt.plot(x, p_numerical, colors[1]+l)

plt.xlabel(r'$L$')
plt.ylabel(r'$p$')
plt.xlim(x[0], x[-1])
plt.legend()
plt.grid()
plt.savefig(join('images', 'pipes', 'nonstationary.pdf'), transparent=True)
