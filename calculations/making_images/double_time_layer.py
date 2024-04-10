import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'double_time_layer'
data = join('data', name)
images = join('images', 'diplom', name)

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif sigma
M = 400
tau = 0.01
ss = [0.5, 1.0, 1.5]

d = []
for s in ss:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{1.0}_sigma{s}.npz')))

for f in ['h', 'u']:
    plotting.dif_times(d[0]['x'], d[0][f], d[0]['x_e'], d[0][f'{f}_e'], fr'${f}$', [fr'$t = {t}$' for t in d[0]['ts']],
                       join(images, f'{name}_{f}_s{ss[0]}.png'))
    plotting.dif_times_dif_param(d[0]['x'], (d[1][f], d[2][f]), d[0]['x_e'], d[0][f'{f}_e'], fr'${f}$', [fr'$t = {t}$' for t in d[0]['ts']],
                    join(images, f'{name}_{f}_s{ss[1]}_s{ss[2]}.png'))
    
    # plotting.dif_times_dif_s(d[0]['x'], (d[0][f], d[1][f], d[2][f]), d[0]['x_e'], d[0][f'{f}_e'], fr'${f}$', [fr'$t = {t}$' for t in d[0]['ts']],
    #                 join(images, f'{name}_{f}.png'))

# dif M
Ms = [200, 400, 800]
tau = 0.01
ss = [0.5, 1.0, 1.5]

d = []
for s in ss:
    for M in Ms:
        d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{1.0}_sigma{s}.npz')))

t = d[0]['t']
for f in ['h', 'u']:
    yss = [[d[i+j][f'err_{f}'] for i in range(3)] for j in range(0, 9, 3)]
    ylabel = fr'$\varepsilon_{f}$'
    names = [fr'$M = {M}$' for M in Ms]
    fname = join(images, f'{name}_err_M_{f}.png')

    plotting.err_dif_M(t, yss, ylabel, names, fname)

# dif tau
M = 400
taus = [0.005, 0.01, 0.02]
ss = [0.5, 1.0, 1.5]

d = []
for s in ss:
    for tau in taus:
        d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{1.0}_sigma{s}.npz')))

for f in ['h', 'u']:
    ts = [d[i]['t'] for i in range(3)]
    yss = [[d[i+j][f'err_{f}'] for i in range(3)] for j in range(0, 9, 3)]
    ylabel = fr'$\varepsilon_{f}$'
    names = [fr'$\tau = {tau}$' for tau in taus]
    fname = join(images, f'{name}_err_t_{f}.png')

    plotting.err_dif_tau(ts, yss, ylabel, names, fname)

