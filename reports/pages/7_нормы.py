import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

norms = join(dirname(dirname(__file__)), 'images', 'norms')

r'''
## Нормы
'''

with st.columns(3)[1]:
    mesh_size = st.select_slider(r'$m$', (100, 200, 400, 800), 200)
    tau = st.select_slider(r'$\tau$', (0.01, 0.005, 0.0025), 0.005)
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0)

cols = st.columns(2)
with cols[0]:
    fname = join(norms, f'norm_h_ms{mesh_size}_tau{tau}_s{s}.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(norms, f'errnorm_h_ms{mesh_size}_tau{tau}_s{s}.png')
    if exists(fname):
            st.image(Image.open(fname))

cols = st.columns(2)
with cols[0]:
    fname = join(norms, f'norm_u_ms{mesh_size}_tau{tau}_s{s}.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(norms, f'errnorm_u_ms{mesh_size}_tau{tau}_s{s}.png')
    if exists(fname):
            st.image(Image.open(fname))
