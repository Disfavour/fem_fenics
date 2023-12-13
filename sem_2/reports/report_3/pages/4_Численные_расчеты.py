import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

images = join(dirname(dirname(__file__)), 'images')

r'''
## Численные расчеты

### Решение на различные моменты времени
'''
with st.columns(3)[1]:
    hu = st.selectbox('Решение', ('h', 'u'))
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1, 1.5), 1, key=1)
    tau = st.select_slider(r'$\tau$', (0.0025, 0.005, 0.01), 0.005, key=2)
    m = st.select_slider(r'$m$', (100, 200, 400), 200, key=3)

cols = st.columns(2)
with cols[0]:
    r'''#### Тестовая задача 1'''
    fname = join(images, f't1_{hu}_dif_t_s{s}_tau{tau}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))
with cols[1]:
    r'''#### Тестовая задача 2'''
    fname = join(images, f't2_{hu}_dif_t_s{s}_tau{tau}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))

r'''
### Полная механическая энергия при разных сетках
'''
with st.columns(3)[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1, 1.5), 1, key=4)
    tau = st.select_slider(r'$\tau$', (0.0025, 0.005, 0.01), 0.005, key=5)

cols = st.columns(2)
with cols[0]:
    r'''#### Тестовая задача 1'''
    fname = join(images, f't1_E_dif_m_s{s}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
with cols[1]:
    r'''#### Тестовая задача 2'''
    fname = join(images, f't2_E_dif_m_s{s}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))

r'''
### Полная механическая энергия при разных $\sigma$
'''
with st.columns(3)[1]:
    tau = st.select_slider(r'$\tau$', (0.0025, 0.005, 0.01), 0.005, key=6)
    m = st.select_slider(r'$m$', (100, 200, 400), 200, key=7)

cols = st.columns(2)
with cols[0]:
    r'''#### Тестовая задача 1'''
    fname = join(images, f't1_E_dif_s_tau{tau}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))
with cols[1]:
    r'''#### Тестовая задача 2'''
    fname = join(images, f't2_E_dif_s_tau{tau}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))

r'''
### Полная механическая энергия при разных $\tau$
'''
with st.columns(3)[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1, 1.5), 1, key=8)
    m = st.select_slider(r'$m$', (100, 200, 400), 200, key=9)

cols = st.columns(2)
with cols[0]:
    r'''#### Тестовая задача 1'''
    fname = join(images, f't1_E_dif_tau_s{s}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))
with cols[1]:
    r'''#### Тестовая задача 2'''
    fname = join(images, f't2_E_dif_tau_s{s}_m{m}.png')
    if exists(fname):
        st.image(Image.open(fname))
