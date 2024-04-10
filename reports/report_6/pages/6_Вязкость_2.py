import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name_ = 'viscosity2'
mesh_size = 200
tau = 0.0025
tests = list(range(1, 4))
degrees = list(range(1, 4))
kappas = [0] + [0.01 * 2**i for i in range(8)]

dir = join(dirname(dirname(dirname(dirname(__file__)))), 'images', 'monotonization', 'linear', name_)

'## Вязкость 2'

name = st.sidebar.radio("Содержание", (
    "Вычислительная схема",
    "Численные эксперименты",
    ))

if name == "Вычислительная схема":
    r'''
Монотонизирующая схема
$$
\left( 1 + \kappa \tau A \right) \frac {u^{n+1} - u^n} \tau + A \frac {u^{n+1} + u^n} 2 = 0
$$

Схема с вязкостью 1
$$
\frac {u^{n+1} - u^n} \tau + (A + \kappa \tau A^* A) \frac {u^{n+1} + u^n} 2 = 0
$$

Приведем добавок монотонизирующей схемы к виду вязкости из схемы 1
$$
(1 + \kappa \tau A^* A) \frac {u^{n+1} - u^n} \tau + A \frac {u^{n+1} + u^n} 2 = 0
$$
'''

elif name == "Численные эксперименты":
    "### Численные эксперименты"
    with st.columns([0.3, 0.4, 0.3])[1]:
        test = st.selectbox('Тестовая задача', ('1', '2', '3'))
        cols = st.columns(2)
        with cols[0]:
            k = st.select_slider(r'$\kappa$', kappas)
        with cols[1]:
            p = st.select_slider(r'$p$', degrees)

    cols = st.columns(2)
    with cols[0]:
        '#### Решения на различные моменты времени'
        fname = join(dir, f'k{k}_p{p}_test{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Закон сохранения массы'
        fname = join(dir, f'm_k{k}_p{p}_test{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    