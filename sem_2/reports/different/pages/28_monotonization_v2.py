import streamlit as st
from PIL import Image
from os.path import dirname, join, exists
import numpy as np

name = 'monotonization_5'
mesh_sizes = (400,)
taus = (0.01,)
alfas = [i * 500 for i in range(0, 21)]
gammas = np.arange(0, 5.1, 0.25)
mesh_sizes = (400,)
taus = (0.01,)

mesh_size = mesh_sizes[0]
tau = taus[0]

hl, hr = 2, 1
dir = join(dirname(dirname(__file__)), 'images', name)

rf'## Результаты решения начально-краевой задачи {name}'
'### Решения на различные моменты времени'
cols = st.columns(2)
with cols[0]:
    alfa1 = st.select_slider(r'$alfa$', alfas, key=1)
    gamma1 = st.select_slider(r'$gamma$', gammas, key=2)
with cols[1]:
    alfa2 = st.select_slider(r'$alfa$', alfas, key=3)
    gamma2 = st.select_slider(r'$gamma$', gammas, key=4)
cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa1}_gamma{gamma1}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa2}_gamma{gamma2}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_u_ms{mesh_size}_tau{tau}_alfa{alfa1}_gamma{gamma1}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_u_ms{mesh_size}_tau{tau}_alfa{alfa2}_gamma{gamma2}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_alfa{alfa1}_gamma{gamma1}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_alfa{alfa2}_gamma{gamma2}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_alfa{alfa1}_gamma{gamma1}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_alfa{alfa2}_gamma{gamma2}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')


