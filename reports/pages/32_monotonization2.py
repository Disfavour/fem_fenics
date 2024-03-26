import streamlit as st
from PIL import Image
from os.path import dirname, join, exists
import numpy as np

name = 'monotonization_2'
mesh_sizes = (800, 400, 200)
taus = (0.005, 0.01, 0.02)
alfas = [-(10 ** i) for i in range(4, 0, -1)] + list(range(-9, -1)) + [round(i*0.1 - 1, 2) for i in range(21)] + list(range(2, 10)) + [10 ** i for i in range(1, 5)]
gammas = list(np.arange(-5, 5.1, 0.5))

hl, hr = 2, 1
dir = join(dirname(dirname(__file__)), 'images', name)

rf'## Результаты решения начально-краевой задачи {name}'
'### Решения на различные моменты времени'

mesh_size = st.select_slider(r'$M$', reversed(mesh_sizes))
tau = st.select_slider(r'$\tau$', taus)

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
