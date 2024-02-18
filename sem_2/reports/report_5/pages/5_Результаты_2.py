import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

dir = join(dirname(dirname(__file__)), 'images')
hl, hr = 2, 1

rf'## Результаты решения начально-краевой задачи с параметром $h_l={hl}$'

'### Решения на различные моменты времени'
cols = st.columns(3)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=1)
with cols[1]:
    tau = st.select_slider(r'$\tau$', (0.00125, 0.0025, 0.005, 0.01), 0.005, key=2)
with cols[2]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=3)
cols = st.columns(2)
with cols[0]:
    '#### Высота'
    fname = join(dir, f'hl{hl}_hr{hr}_h_ms{ms}_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    '#### Скорость'
    fname = join(dir, f'hl{hl}_hr{hr}_u_ms{ms}_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

'### Влияние сеток по пространству'
cols = st.columns(2)
with cols[0]:
    tau = st.select_slider(r'$\tau$', (0.00125, 0.0025, 0.005, 0.01), 0.005, key=4)
with cols[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=5)
cols = st.columns(2)
with cols[0]:
    '#### Норма ошибки высоты'
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    '#### Норма ошибки скорости'
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with st.columns([0.25, 0.5, 0.25])[1]:
    '#### Полная механическая энергия'
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

'### Влияние шага по времени'
cols = st.columns(2)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=6)
with cols[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=7)
cols = st.columns(2)    
with cols[0]:
    '#### Норма ошибки высоты'
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    '#### Норма ошибки скорости'
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with st.columns([0.25, 0.5, 0.25])[1]:
    '#### Полная механическая энергия'
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

'### Влияние весового параметра'
cols = st.columns(2)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=8)
with cols[1]:
    tau = st.select_slider(r'$\tau$', (0.00125, 0.0025, 0.005, 0.01), 0.005, key=9)
cols = st.columns(2)
with cols[0]:
    '#### Норма ошибки высоты'
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    '#### Норма ошибки скорости'
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with st.columns([0.25, 0.5, 0.25])[1]:
    '#### Полная механическая энергия'
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
