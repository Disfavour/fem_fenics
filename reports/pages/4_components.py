import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

images = join(dirname(dirname(__file__)), 'images')

r'''
## Компоненты
'''

with st.columns(3)[1]:
    m = st.select_slider(r'$m$', (100, 200, 400), 200)
    tau = st.select_slider(r'$\tau$', (0.0025, 0.005, 0.01), 0.005)
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1, 1.5), 1)

# ms100_tau0.01_s0.75_c0.png
cols = st.columns(2)
with cols[0]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c0.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c0_dif.png')
    if exists(fname):
            st.image(Image.open(fname))

cols = st.columns(2)
with cols[0]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c1.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c1_dif.png')
    if exists(fname):
            st.image(Image.open(fname))

cols = st.columns(2)
with cols[0]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c2.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c2_dif.png')
    if exists(fname):
            st.image(Image.open(fname))
    
cols = st.columns(2)
with cols[0]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c3.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c3_dif.png')
    if exists(fname):
            st.image(Image.open(fname))

cols = st.columns(2)
with cols[0]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c4.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'ms{m}_tau{tau}_s{s}_c4_dif.png')
    if exists(fname):
            st.image(Image.open(fname))

# ms{ms}_tau{tau}_s{s}_m.png
fname = join(images, f'ms{m}_tau{tau}_s{s}_m.png')
if exists(fname):
        st.image(Image.open(fname))

fname = join(images, f'ms{m}_tau{tau}_s{s}_E.png')
if exists(fname):
        st.image(Image.open(fname))
