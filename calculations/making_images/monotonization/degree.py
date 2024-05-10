import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'double_time_layer'
data = join('data', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif sigma
M = 200
tau = 0.005
s = 0.6
degree = [2, 1]

d = []
for p in degree:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{p}_sigma{s}.npz')))

for f in ['h', 'u']:
    xs = [d[0]['x_e']] + [i['x'] for i in d]
    yss = []
    for t in range(d[0]['ts'].size):
        ys = [d[0][f'{f}_e'][t]]
        for i in d:
            ys.append(i[f][t])
        yss.append(ys)
        
    ylabel = fr'${f}$'
    labels = [fr'$t = {t}$' for t in d[0]['ts']]
    fname = join(images, f'monotonization_p_{f}.pdf')
    plotting.dif_times(xs, yss, ylabel, labels, fname)

# dif
Ms = [200]
tau = 0.005
degree = [1, 2]

d = []

for p in degree:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{p}_sigma{s}.npz')))

for f in ['h', 'u']:
    ts = [i['t'] for i in d]
    yss = [[i[f'err_{f}']] for i in d]
    ylabel = fr'$\varepsilon_{f}$'
    labels = [fr'$p = {p}$' for p in degree]
    fname = join(images, f'monotonization_p_err_{f}.pdf')
    plotting.dif_params(ts, yss, ylabel, labels, fname)
