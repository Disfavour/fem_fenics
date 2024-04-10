import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name_ = 'triple_time_layer_1'
mesh_sizes = (800, 400, 200)
taus = (0.005, 0.01, 0.02)
thetas = (0.5, 1.0, 1.25, 1.5)
sigmas = (1.0,)

sigma = sigmas[0]
mesh_sizes = tuple(reversed(mesh_sizes))

hl, hr = 2, 1
dir = join(dirname(dirname(dirname(dirname(__file__)))), 'images', name_)

r'''
## Трехслойная схема с весами 1
'''

name = st.sidebar.radio("Содержание", (
    "Вычислительная схема",
    "Численные решения",
    'Влияние сеток по пространству',
    'Влияние шага по времени',
    'Влияние весового параметра $\sigma$'
    ))

if name == "Вычислительная схема":
    r'''
Аппроксимация по времени
$$
\frac {\partial \varphi} {\partial t} \approx \theta \frac {\varphi^{n+1} - \varphi^n} {\tau} + (1 - \theta) \frac {\varphi^{n} - \varphi^{n-1}} {\tau},
\quad \theta = \operatorname{const}
$$

Неявная схема с весами
$$
\begin{aligned}
    & \theta \frac {h^{n+1} - h^n} {\tau} + (1 - \theta) \frac {h^{n} - h^{n-1}} {\tau} + \frac {\partial} {\partial x} (h u)^{n+1} = 0
    \\[0.5 cm]
    & \theta \frac {(h u)^{n+1} - (h u)^n} {\tau} + (1 - \theta) \frac {(h u)^{n} - (h u)^{n-1}} {\tau}
    + \frac {\partial} {\partial x} (h u^2)^{n+1} + g h^{n+1} \frac {\partial h^{n+1}} {\partial x} = 0
\end{aligned}
$$

- $\theta = \displaystyle \frac 1 2$ — предельное по устойчивости
- $\theta = 1$ — двухслойная чисто неявная схема
- $\theta = \displaystyle \frac 3 2$ — схема второго порядка аппроксимации по $\tau$
    '''

elif name == "Численные решения":
    '### Решения на различные моменты времени'
    cols = st.columns(3)
    with cols[0]:
        mesh_size = st.select_slider(r'$M$', mesh_sizes, key=1)
    with cols[1]:
        tau = st.select_slider(r'$\tau$', taus, key=2)
    with cols[2]:
        theta = st.select_slider(r'$\theta$', thetas, key=3)

    cols = st.columns(2)
    with cols[0]:
        '#### Высота'
        print(join(dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.png'))
        fname = join(dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Скорость'
        fname = join(dir, f'hl{hl}_hr{hr}_u_ms{mesh_size}_tau{tau}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')

elif name == 'Влияние сеток по пространству':
    '### Влияние сеток по пространству'
    cols = st.columns(2)
    with cols[0]:
        tau = st.select_slider(r'$\tau$', taus, key=4)
    with cols[1]:
        theta = st.select_slider(r'$\theta$', thetas, key=5)
    cols = st.columns(2)
    with cols[0]:
        '#### Норма ошибки высоты'
        fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Норма ошибки скорости'
        fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with st.columns([0.25, 0.5, 0.25])[1]:
        '#### Полная механическая энергия'
        fname = join(dir, f'hl{hl}_hr{hr}_E_dif_ms_tau{tau}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')

elif name == 'Влияние шага по времени':
    '### Влияние шага по времени'
    cols = st.columns(2)
    with cols[0]:
        mesh_size = st.select_slider(r'$M$', mesh_sizes, key=6)
    with cols[1]:
        theta = st.select_slider(r'$\theta$', thetas, key=7)
    cols = st.columns(2)    
    with cols[0]:
        '#### Норма ошибки высоты'
        fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Норма ошибки скорости'
        fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with st.columns([0.25, 0.5, 0.25])[1]:
        '#### Полная механическая энергия'
        fname = join(dir, f'hl{hl}_hr{hr}_E_dif_tau_ms{mesh_size}_theta{theta}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')

elif name == 'Влияние весового параметра $\sigma$':
    r'### Влияние весового параметра $\theta$'
    cols = st.columns(2)
    with cols[0]:
        mesh_size = st.select_slider(r'$M$', mesh_sizes, key=8)
    with cols[1]:
        tau = st.select_slider(r'$\tau$', taus, key=9)
    cols = st.columns(2)
    with cols[0]:
        '#### Норма ошибки высоты'
        fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Норма ошибки скорости'
        fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with st.columns([0.25, 0.5, 0.25])[1]:
        '#### Полная механическая энергия'
        fname = join(dir, f'hl{hl}_hr{hr}_E_dif_theta_ms{mesh_size}_tau{tau}_sigma{sigma}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
