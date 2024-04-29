import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'triple_time_layer_1'
data = join('data', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif theta
M = 400
tau = 0.01
thetas = (0.5, 1.0, 1.5)
s = 1.0

d = []
for theta in thetas:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{theta}_sigma{s}.npz')))

for f in ['h', 'u']:
    xs = [d[0]['x_e']] + [i['x'] for i in d]
    yss = []
    for t in range(d[0]['ts'].size):
        ys = []
        for i in d:
            ys.append(i[f][t])
        ys.reverse()
        ys = [d[0][f'{f}_e'][t]] + ys
        yss.append(ys)
        
    ylabel = fr'${f}$'
    labels = [fr'$t = {t}$' for t in d[0]['ts']]
    fname = join(images, f'{name}_{f}.pdf')
    plotting.dif_times(xs, yss, ylabel, labels, fname)

# dif M
Ms = [200, 400, 800]
tau = 0.01
thetas = (0.5, 1.0, 1.5)
s = 1.0

d = []
for M in Ms:
    tmp = []
    for theta in thetas:
        tmp.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{theta}_sigma{s}.npz')))
    d.append(tmp)

for f in ['h', 'u']:
    ts = [i[0]['t'] for i in d]
    yss = [[j[f'err_{f}'] for j in i] for i in d]
    ylabel = fr'$\varepsilon_{f}$'
    labels = [fr'$M = {M}$' for M in Ms]
    fname = join(images, f'{name}_err_M_{f}.pdf')
    plotting.dif_params(ts, yss, ylabel, labels, fname)

# dif tau
M = 400
taus = [0.005, 0.01, 0.02]
taus.reverse()
thetas = (0.5, 1.0, 1.5)
s = 1.0

d = []
for tau in taus:
    tmp = []
    for theta in thetas:
        tmp.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{theta}_sigma{s}.npz')))
    d.append(tmp)

for f in ['h', 'u']:
    ts = [i[0]['t'] for i in d]
    yss = [[j[f'err_{f}'] for j in i] for i in d]
    ylabel = fr'$\varepsilon_{f}$'
    labels = [fr'$\tau = {tau}$' for tau in taus]
    fname = join(images, f'{name}_err_tau_{f}.pdf')
    plotting.dif_params(ts, yss, ylabel, labels, fname)
