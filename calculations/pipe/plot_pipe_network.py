import numpy as np
import pyvista as pv
from os.path import dirname, join

data = np.load(join('data', 'diamond_network.npy'))
p = data[::2, :]
m = data[1::2, :]


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


l1 = np.row_stack((points[0], points[1], points[3], points[5], points[6]))
l2 = np.row_stack((points[0], points[1], points[4], points[6]))
l3 = np.row_stack((points[0], points[2], points[5], points[6]))

p1 = np.concatenate((p[0], p[1], p[3], p[5], p[6]))
p2 = np.concatenate((p[0], p[1], p[4], p[6]))
p3 = np.concatenate((p[0], p[2], p[5], p[6]))

m1 = np.concatenate((m[0], m[1], m[3], m[5], m[6]))
m2 = np.concatenate((m[0], m[1], m[4], m[6]))
m3 = np.concatenate((m[0], m[2], m[5], m[6]))

ps = [p1, p2, p3]
ms = [m1, m2, m3]

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

for fs, fname in zip([ps, ms], ['p', 'm']):
    pl = pv.Plotter(off_screen=True)
    for points, f in zip([l1, l2, l3], fs):
        spline = pv.Spline(points)
        spline[''] = f
        tube = spline.tube(radius=0.05)
        pl.add_mesh(tube, smooth_shading=True, scalar_bar_args=dict(n_labels=3))

    #bar = pl.add_scalar_bar(n_labels=2)
    #bar.GetTitleTextProperty().SetLineSpacing(1.5)

    #pl.show_axes()
    pl.set_viewup([0, 1, 0])
    pl.camera.zoom(1.8)
    #pl.save_graphic("img.pdf")
    pl.screenshot(fname, transparent_background=True, window_size=(1920, 1080), scale=2)
    #pl.show()
