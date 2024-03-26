import streamlit as st
from PIL import Image
from os.path import dirname, join, exists
import numpy as np

name = 'monotonization_artificial_viscosity'
mesh_sizes = (800, 400, 200)
taus = (0.005, 0.01, 0.02)
degrees = list(range(1, 4))
T = 7.0
sigma = 0.5
alfas = [0] + [2**i for i in range(8)]
gammas = list(range(4))
kappas = [0.5]
ntests = (1, 2)

kappa = kappas[0]

mesh_sizes = list(reversed(mesh_sizes))

dir = join(dirname(dirname(__file__)), 'images', name)

rf'## Результаты решения начально-краевой задачи {name}'
'### Решения на различные моменты времени'

with st.columns([0.25, 0.5, 0.25])[1]:
    ntest = st.select_slider(r'$ntest$', ntests)
    mesh_size = st.select_slider(r'$M$', mesh_sizes)
    tau = st.select_slider(r'$\tau$', taus)

cols = st.columns(2)
with cols[0]:
    degree1 = st.select_slider(r'$degree$', degrees, key=1)
    alfa1 = st.select_slider(r'$alfa$', alfas, key=2)
    gamma1 = st.select_slider(r'$gamma$', gammas, key=3)
with cols[1]:
    degree2 = st.select_slider(r'$degree$', degrees, key=5)
    alfa2 = st.select_slider(r'$alfa$', alfas, key=6)
    gamma2 = st.select_slider(r'$gamma$', gammas, key=7)

cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'ms{mesh_size}_tau{tau}_degree{degree1}_alfa{alfa1}_gamma{gamma1}_kappa{kappa}_ntest{ntest}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'ms{mesh_size}_tau{tau}_degree{degree2}_alfa{alfa2}_gamma{gamma2}_kappa{kappa}_ntest{ntest}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
