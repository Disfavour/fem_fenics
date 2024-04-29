import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


images = join('images', 'diplom')
for name in ['weight', 'viscosity1', 'viscosity2']:
    data = join('data', 'monotonization', 'linear', name)

    makedirs(data, exist_ok=True)
    makedirs(images, exist_ok=True)

    # dif kappa
    tests = list(range(1, 4))
    ks = [0.0, 0.04, 0.08]
    for test in tests:
        d = []
        for k in ks:
            d.append(np.load(join(data, f'k{k}_test{test}.npz')))

        f = 'u'
        x = d[0]['x']
        yss = [i[f] for i in d]
            
        ylabel = fr'${f}$'
        labels = [fr'$\kappa = {k}$' for k in ks]
        fname = join(images, f'{name}_test{test}.pdf')
        plotting.dif_times_special(x, yss, ylabel, labels, fname)

        ts = [i['t'] for i in d]
        ys = [i['err'] for i in d]
        xlabel = r'$t$'
        ylabel = fr'$\varepsilon_{f}$'
        fname = join(images, f'{name}_err_test{test}.pdf')
        plotting.base(ts, ys, xlabel, ylabel, labels, fname)
