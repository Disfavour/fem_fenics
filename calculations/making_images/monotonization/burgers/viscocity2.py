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

test = 2
ks = [0.02, 0.04, 0.08]

d = []
for k in ks:
    d.append(np.load(join(data, f'k{k}_test{test}.npz')))

xs = [d[0]['x_e']] + [i['x'] for i in d]
yss = []
for t in range(d[0]['ts'].size):
    ys = []
    for i in d:
        ys.append(i['u'][t])
    ys.reverse()
    ys = [d[0][f'u_e'][t]] + ys
    yss.append(ys)
    
ylabel = fr'$u$'
labels = [fr'$t = {t}$' for t in d[0]['ts']]
fname = join(images, f'burgers_{name}.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)
