import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'nonlinear_viscosity'
data = join('data', 'monotonization', 'linear', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

test = 1

M = 200
tau = 0.0025

a = 1
gammas = [0, 1, 2]

d = []
for y in gammas:
    d.append(np.load(join(data, f'a{a}_y{y}_test{test}.npz')))

xs = [i['x'] for i in d]
yss = []
for t in range(d[0]['ts'].size):
    ys = []
    for i in d:
        ys.append(i['u'][t])
    ys.reverse()
    yss.append(ys)
    
ylabel = fr'$u$'
labels = [fr'$t = {t}$' for t in d[0]['ts']]
fname = join(images, f'transport_{name}_dif_y.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)


# dif M
Ms = [100, 200, 400]
tau = 0.0025

d = []
for M in Ms:
    tmp = []
    for y in gammas:
        tmp.append(np.load(join(data, f'M{M}_tau{tau}_a{a}_y{y}.npz')))
    d.append(tmp)

f = 'u'
ts = [i[0]['t'] for i in d]
yss = [[j[f'err'] for j in i] for i in d]
ylabel = fr'$\varepsilon_{f}$'
labels = [fr'$M = {M}$' for M in Ms]
fname = join(images, f'transport_{name}_dif_y_err_M.pdf')
plotting.dif_params(ts, yss, ylabel, labels, fname)

# dif tau
M = 200
taus = [0.005, 0.0025, 0.00125]

d = []
for tau in taus:
    tmp = []
    for y in gammas:
        tmp.append(np.load(join(data, f'M{M}_tau{tau}_a{a}_y{y}.npz')))
    d.append(tmp)

f = 'u'
ts = [i[0]['t'] for i in d]
yss = [[j[f'err'] for j in i] for i in d]
ylabel = fr'$\varepsilon_{f}$'
labels = [fr'$\tau = {tau}$' for tau in taus]
fname = join(images, f'transport_{name}_dif_y_err_tau.pdf')
plotting.dif_params(ts, yss, ylabel, labels, fname)
