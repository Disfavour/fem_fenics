import numpy as np
import pyvista as pv
from os.path import dirname, join

data = np.load(join('data', 'diamond_network_components.npy'))
p = data[::2, :]
m = data[1::2, :]

p_nodes = [p[0][0], p[0][-1], p[1][-1], p[5][0], p[6][0], p[6][-1]]
m_nodes = [m[0][0], m[0][-1], m[1][-1], m[5][0], m[6][0], m[6][-1]]

p_nodes = np.array((p_nodes, p_nodes, p_nodes)).T
m_nodes = np.array((m_nodes, m_nodes, m_nodes)).T

n = data.shape[1]   # 51
L = 1

points = []
for i in range(7):
    points.append(np.zeros((n, 3)))

points[0][:, 0] = np.linspace(0, L, n)

points[1][:, 0] = np.linspace(L, L*(1 + np.sqrt(3)/2), n)
points[1][:, 1] = np.linspace(0, L/2, n)

points[2][:, 0] = np.linspace(L, L*(1 + np.sqrt(3)/2), n)
points[2][:, 1] = np.linspace(0, -L/2, n)

points[3][:, 0] = L*(1 + np.sqrt(3)/2)
points[3][:, 1] = np.linspace(L/2, -L/2, n)

points[4][:, 0] = np.linspace(L*(1 + np.sqrt(3)/2), L*(1 + np.sqrt(3)), n)
points[4][:, 1] = np.linspace(L/2, 0, n)

points[5][:, 0] = np.linspace(L*(1 + np.sqrt(3)/2), L*(1 + np.sqrt(3)), n)
points[5][:, 1] = np.linspace(-L/2, 0, n)

points[6][:, 0] = np.linspace(L*(1 + np.sqrt(3)), L*(2 + np.sqrt(3)), n)

r = 0.05
z0 = r - 0.001
# узлы
nodes = []
nodes.append(np.array(( # 0
    (0, 0, r),
    (0, 0, z0)
)))
nodes.append(np.array(( # 1
    (L, 0, r),
    (L, 0, z0)
)))
nodes.append(np.array(( # 2
    (L*(1 + np.sqrt(3)/2), L/2, r),
    (L*(1 + np.sqrt(3)/2), L/2, z0)
)))
nodes.append(np.array(( # 3
    (L*(1 + np.sqrt(3)/2), -L/2, r),
    (L*(1 + np.sqrt(3)/2), -L/2, z0)
)))
nodes.append(np.array(( # 4
    (L*(1 + np.sqrt(3)), 0, r),
    (L*(1 + np.sqrt(3)), 0, z0)
)))
nodes.append(np.array(( # 5
    (L*(2 + np.sqrt(3)), 0, r),
    (L*(2 + np.sqrt(3)), 0, z0)
)))

pv.global_theme.font.family = 'arial'
# pv.global_theme.font.size = 20
# pv.global_theme.font.title_size = 40
pv.global_theme.font.label_size = 54
pv.global_theme.font.color = 'black'
pv.global_theme.font.fmt = '%.4g'

pv.global_theme.colorbar_horizontal.position_x = 0.1
pv.global_theme.colorbar_horizontal.position_y = 0.05
pv.global_theme.colorbar_horizontal.width = 0.8
pv.global_theme.colorbar_horizontal.height = 0.1

for fs, f_nodes, fname in zip([p, m], [p_nodes, m_nodes], ['p_components', 'm_components']):
    pl = pv.Plotter(off_screen=True)
    for pipe, f in zip(points, fs):
        spline = pv.Spline(pipe)
        spline[''] = f
        tube = spline.tube(radius=r)
        pl.add_mesh(tube, smooth_shading=True, scalar_bar_args=dict(n_labels=3))

    for node, f in zip(nodes, f_nodes):
        spline = pv.Spline(node)
        spline[''] = f
        tube = spline.tube(radius=2*r)
        pl.add_mesh(tube, smooth_shading=True, scalar_bar_args=dict(n_labels=3))

    pl.set_viewup([0, 1, 0])
    pl.camera.zoom(1.8)
    pl.screenshot(fname, transparent_background=True, window_size=(1920, 1080), scale=2)
