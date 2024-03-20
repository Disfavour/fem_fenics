import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name = 'triple_time_layer_1'
mesh_sizes = (800, 400, 200)
taus = (0.005, 0.01, 0.02)
thetas = (0.5, 1.0, 1.25, 1.5)
sigmas = (1.0,)

sigma = sigmas[0]
mesh_sizes = tuple(reversed(mesh_sizes))

hl, hr = 2, 1
dir = join(dirname(dirname(__file__)), 'images', name)

rf'## Результаты решения начально-краевой задачи {name}'

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
