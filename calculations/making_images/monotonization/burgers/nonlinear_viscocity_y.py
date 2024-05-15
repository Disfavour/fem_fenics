import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'nonlinear_viscosity'
data = join('data', 'monotonization', 'nonlinear', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# [0, 50, 100, 200, 500, 1000]
a = 1
y_s = [1, 2, 3]

M = 200
tau = 0.0025

d = []
for y in y_s:
    d.append(np.load(join(data, f'M{M}_tau{tau}_a{a}_y{y}.npz')))

xs = [d[0]['x_e']] + [i['x'] for i in d]
yss = []
for t in range(d[0]['ts'].size):
    ys = [d[0][f'u_e'][t]]
    for i in d:
        ys.append(i['u'][t])
    yss.append(ys)
    
ylabel = fr'$u$'
labels = [fr'$t = {t}$' for t in d[0]['ts']]
fname = join(images, f'burgers_{name}_dif_y.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)


y_s = [0] + y_s

# dif M
Ms = [100, 200, 400]
tau = 0.0025

d = []
for M in Ms:
    tmp = []
    for y in y_s:
        tmp.append(np.load(join(data, f'M{M}_tau{tau}_a{a}_y{y}.npz')))
    d.append(tmp)

f = 'u'
ts = [i[0]['t'] for i in d]
yss = [[j[f'err'] for j in i] for i in d]
ylabel = fr'$\varepsilon_{f}$'
labels = [fr'$M = {M}$' for M in Ms]
fname = join(images, f'burgers_{name}_dif_y_err_M.pdf')
plotting.dif_params(ts, yss, ylabel, labels, fname)

# dif tau
M = 200
taus = [0.005, 0.0025, 0.00125]

d = []
for tau in taus:
    tmp = []
    for y in y_s:
        tmp.append(np.load(join(data, f'M{M}_tau{tau}_a{a}_y{y}.npz')))
    d.append(tmp)

f = 'u'
ts = [i[0]['t'] for i in d]
yss = [[j[f'err'] for j in i] for i in d]
ylabel = fr'$\varepsilon_{f}$'
labels = [fr'$\tau = {tau}$' for tau in taus]
fname = join(images, f'burgers_{name}_dif_y_err_tau.pdf')
plotting.dif_params(ts, yss, ylabel, labels, fname)
