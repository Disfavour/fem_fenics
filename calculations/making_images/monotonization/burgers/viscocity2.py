import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'viscosity2'
data = join('data', 'monotonization', 'nonlinear', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

M = 200
tau = 0.0025

ks = [0.02, 0.04, 0.08]

d = []
for k in ks:
    d.append(np.load(join(data, f'M{M}_tau{tau}_k{k}.npz')))

xs = [d[0]['x_e']] + [i['x'] for i in d]
yss = []

for t in range(d[0]['ts'].size):
    ys = [d[0][f'u_e'][t]]
    for i in d:
        ys.append(i['u'][t])
    yss.append(ys)
    
ylabel = fr'$u$'
labels = [fr'$t = {t}$' for t in d[0]['ts']]
fname = join(images, f'burgers_{name}.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)

ks = [0.0] + ks

# # dif M
# Ms = [100, 200, 400]
# tau = 0.0025

# d = []
# for M in Ms:
#     tmp = []
#     for k in ks:
#         tmp.append(np.load(join(data, f'M{M}_tau{tau}_k{k}.npz')))
#     d.append(tmp)

# f = 'u'
# ts = [i[0]['t'] for i in d]
# yss = [[j[f'err'] for j in i] for i in d]
# ylabel = fr'$\varepsilon_{f}$'
# labels = [fr'$M = {M}$' for M in Ms]
# fname = join(images, f'burgers_{name}_err_M.pdf')
# plotting.dif_params(ts, yss, ylabel, labels, fname)

# # dif tau
# M = 200
# taus = [0.005, 0.0025, 0.00125]

# d = []
# for tau in taus:
#     tmp = []
#     for k in ks:
#         tmp.append(np.load(join(data, f'M{M}_tau{tau}_k{k}.npz')))
#     d.append(tmp)

# f = 'u'
# ts = [i[0]['t'] for i in d]
# yss = [[j[f'err'] for j in i] for i in d]
# ylabel = fr'$\varepsilon_{f}$'
# labels = [fr'$\tau = {tau}$' for tau in taus]
# fname = join(images, f'burgers_{name}_err_tau.pdf')
# plotting.dif_params(ts, yss, ylabel, labels, fname)
