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
alfas = [0, 16, 32]
y = 0

d = []
for a in alfas:
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
fname = join(images, f'transport_{name}_dif_a.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)

#
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
