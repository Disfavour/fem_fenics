import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


images = join('images', 'diplom')
name = 'nonlinear_viscosity'
data = join('data', 'monotonization', 'linear', name)

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif kappa
tests = list(range(1, 4))
alfas = [0, 16, 32]
y = 0#list(range(3))
for test in tests:
    d = []
    for a in alfas:
        d.append(np.load(join(data, f'a{a}_y{y}_test{test}.npz')))

    f = 'u'
    x = d[0]['x']
    yss = [i[f] for i in d]
        
    ylabel = fr'${f}$'
    labels = [fr'$\alpha = {a}$' for a in alfas]
    fname = join(images, f'{name}_y{y}_test{test}.pdf')
    plotting.dif_times_special(x, yss, ylabel, labels, fname)

    ts = [i['t'] for i in d]
    ys = [i['err'] for i in d]
    xlabel = r'$t$'
    ylabel = fr'$\varepsilon_{f}$'
    fname = join(images, f'{name}_err_y{y}_test{test}.pdf')
    plotting.base(ts, ys, xlabel, ylabel, labels, fname)

a = 1
gammas = list(range(3))
for test in tests:
    d = []
    for y in gammas:
        d.append(np.load(join(data, f'a{a}_y{y}_test{test}.npz')))

    f = 'u'
    x = d[0]['x']
    yss = [i[f] for i in d]
        
    ylabel = fr'${f}$'
    labels = [fr'$\gamma = {y}$' for y in gammas]
    fname = join(images, f'{name}_a{a}_test{test}.pdf')
    plotting.dif_times_special(x, yss, ylabel, labels, fname)

    ts = [i['t'] for i in d]
    ys = [i['err'] for i in d]
    xlabel = r'$t$'
    ylabel = fr'$\varepsilon_{f}$'
    fname = join(images, f'{name}_err_a{a}_test{test}.pdf')
    plotting.base(ts, ys, xlabel, ylabel, labels, fname)