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

test = 2
alfas = [128, 512, 1024]
y = 0

d = []
for a in alfas:
    d.append(np.load(join(data, f'a{a}_y{y}_test{test}.npz')))

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
fname = join(images, f'burgers_{name}_dif_a.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)

#
a = 1
gammas = [1, 2, 3]

d = []
for y in gammas:
    d.append(np.load(join(data, f'a{a}_y{y}_test{test}.npz')))

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
fname = join(images, f'burgers_{name}_dif_y.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)

# base
a = 0
y = 0

d = [np.load(join(data, f'a{a}_y{y}_test{test}.npz'))]

xs = [d[0]['x_e']] + [i['x'] for i in d]
yss = []
for t in range(d[0]['ts'].size):
    ys = [d[0][f'u_e'][t]]
    for i in d:
        ys.append(i['u'][t])
    yss.append(ys)
    
ylabel = fr'$u$'
labels = [fr'$t = {t}$' for t in d[0]['ts']]
fname = join(images, f'burgers.pdf')
plotting.dif_times(xs, yss, ylabel, labels, fname)