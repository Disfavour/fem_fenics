import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'triple_time_layer_2'
data = join('data', name)
images = join('images', 'diplom')

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif sigma
M = 400
tau = 0.01
theta = 0.5
ss = (0.25, 0.5) # 1.0

d = []
for s in ss:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{theta}_sigma{s}.npz')))

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
    fname = join(images, f'{name}_{f}.pdf')
    plotting.dif_times(xs, yss, ylabel, labels, fname)
