import streamlit as st
from PIL import Image
from os.path import dirname, join, exists
import numpy as np

name1 = 'monotonization_1'
name2 = 'monotonization_2'
name3 = 'monotonization_3'
mesh_sizes = (400,)
taus = (0.01,)
alfas = list(np.linspace(0.1, 1, 10)) + [2, 5] + [10**i for i in range(1, 5)]
gammas = [-i for i in reversed(alfas)] + [0] + alfas

mesh_size = mesh_sizes[0]
tau = taus[0]

hl, hr = 2, 1
dir = join(dirname(dirname(__file__)), 'images')

rf'## Результаты решения начально-краевой задачи {name1}'

'### Решения на различные моменты времени'
cols = st.columns(2)
with cols[0]:
    alfa = st.select_slider(r'$alfa$', alfas, key=1)
with cols[1]:
    gamma = st.select_slider(r'$gamma$', gammas, key=2)
cols = st.columns(2)
with cols[0]:
    'monotonization_2'
    fname = join(dir, name2, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    'monotonization_3'
    fname = join(dir, name3, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

with st.columns([0.25, 0.5, 0.25])[1]:
    'monotonization_1'
    fname = join(dir, name1, f'hl{hl}_hr{hr}_h_ms{mesh_size}_tau{tau}_alfa{alfa}_gamma{gamma}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
