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
degree = [1, 2]

d = []
for p in degree:
    d.append(np.load(join(data, f'hl{2}_hr{1}_ms{M}_tau{tau}_theta{p}_sigma{s}.npz')))

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
    fname = join(images, f'monotonization_p_{f}.pdf')
    plotting.dif_times(xs, yss, ylabel, labels, fname)