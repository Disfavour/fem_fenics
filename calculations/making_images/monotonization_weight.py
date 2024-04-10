import sys
sys.path.append('calculations')
import plotting
from os.path import dirname, join, exists
from os import makedirs
import numpy as np


name = 'weight'
data = join('data', 'monotonization', 'linear', name)
images = join('images', 'diplom', 'monotonization', 'linear', name)

makedirs(data, exist_ok=True)
makedirs(images, exist_ok=True)

# dif kappa

# mesh_size = 200
# tau = 0.0025
tests = list(range(1, 4))
# degrees = list(range(1, 4))
# kappas = [0] + [0.01 * 2**i for i in range(8)]

p = 1
ks = [0, 0.04, 0.16]
test = 1

d = []
for k in ks:
    d.append(np.load(join(data, f'k{k}_p{p}_test{test}.npz')))

x = d[0]['x']
yss = [d[i]['u'] for i in range(3)]
ylabel = r'$u$'
names = [fr'$\kappa = {k}$' for k in ks]
fname = join(images, f'{name}.png')
plotting.dif_times_dif_param_linear(x, yss, names, fname)