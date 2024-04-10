import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name_ = 'triple_time_layer_2'
mesh_sizes = (800, 400, 200)
taus = (0.005, 0.01, 0.02)
thetas = (0.5,)
sigmas = (0.25, 0.5, 0.75, 1.0)

theta = thetas[0]
mesh_sizes = tuple(reversed(mesh_sizes))

hl, hr = 2, 1
dir = join(dirname(dirname(dirname(dirname(__file__)))), 'images', name_)

r'''
## Трехслойная схема с весами 2
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
Аппроксимации по времени
$$
\begin{aligned}
	\frac {\partial \varphi} {\partial t} \approx \theta \frac {\varphi^{n+1} - \varphi^n} {\tau} + (1 - \theta) \frac {\varphi^{n} - \varphi^{n-1}} {\tau}
	\\[0.5 cm]
	\varphi^{n+\sigma} = \sigma \varphi^{n+1} + (1 - 2 \sigma) \varphi^n + \sigma \varphi^{n-1}
\end{aligned}
$$

Схема с весами
$$
\begin{aligned}
	& \theta \frac {h^{n+1} - h^n} {\tau} + (1 - \theta) \frac {h^{n} - h^{n-1}} {\tau} + \frac {\partial} {\partial x} (h u)^{n+\sigma} = 0
	\\[0.5 cm]
	& \theta \frac {(h u)^{n+1} - (h u)^n} {\tau} + (1 - \theta) \frac {(h u)^{n} - (h u)^{n-1}} {\tau}
    + \frac {\partial} {\partial x} (h u^2)^{n+\sigma} + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$

Рассматривается при фиксированном $\theta = \displaystyle \frac 1 2$ (2-й порядок аппроксимации по $\tau$)
    '''

elif name == "Численные решения":
    '### Решения на различные моменты времени'
    cols = st.columns(3)
    with cols[0]:
        mesh_size = st.select_slider(r'$M$', mesh_sizes, key=1)
    with cols[1]:
        tau = st.select_slider(r'$\tau$', taus, key=2)
    with cols[2]:
        sigma = st.select_slider(r'$\sigma$', sigmas, key=3)
    cols = st.columns(2)
    with cols[0]:
        '#### Высота'
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
        sigma = st.select_slider(r'$\sigma$', sigmas, key=5)
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
        sigma = st.select_slider(r'$\sigma$', sigmas, key=7)
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
    r'### Влияние весового параметра $\sigma$'
    cols = st.columns(2)
    with cols[0]:
        mesh_size = st.select_slider(r'$M$', mesh_sizes, key=8)
    with cols[1]:
        tau = st.select_slider(r'$\tau$', taus, key=9)
    cols = st.columns(2)
    with cols[0]:
        '#### Норма ошибки высоты'
        fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Норма ошибки скорости'
        fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with st.columns([0.25, 0.5, 0.25])[1]:
        '#### Полная механическая энергия'
        fname = join(dir, f'hl{hl}_hr{hr}_E_dif_sigma_ms{mesh_size}_tau{tau}_theta{theta}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
