import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name_ = 'nonlinear_viscosity'
mesh_size = 200
tau = 0.0025
tests = list(range(1, 4))
degrees = list(range(1, 4))
alfas = [0] + [2**i for i in range(8)]
gammas = list(range(4))

dir = join(dirname(dirname(dirname(dirname(__file__)))), 'images', 'monotonization', 'linear', name_)

'## Нелинейная вязкость'

name = st.sidebar.radio("Содержание", (
    "Вычислительная схема",
    "Численные эксперименты",
    ))

if name == "Вычислительная схема":
    r'''
$$
\begin{aligned}
    &- a \tau^2 \operatorname{div} \left( b \operatorname{grad} u \right)
    \\[0.5 cm]
    & a = \operatorname{const} > 0
    \\[0.5 cm]
    & b = | \operatorname{grad} u |^{\gamma}
\end{aligned}
$$

Для одномерного уравнения переноса
$$
\frac {u^{n+1} - u^n} \tau + \frac \partial {\partial x} \frac {u^{n+1} + u^n} 2
- a \tau^2 \frac \partial {\partial x} \left( b \frac {\partial u} {\partial x} \right) = 0
$$
'''

elif name == "Численные эксперименты":
    "### Численные эксперименты"
    with st.columns([0.3, 0.4, 0.3])[1]:
        test = st.selectbox('Тестовая задача', ('1', '2', '3'))
        cols = st.columns(3)
        with cols[0]:
            a = st.select_slider(r'$a$', alfas)
        with cols[1]:
            y = st.select_slider(r'$\gamma$', gammas)
        with cols[2]:
            p = st.select_slider(r'$p$', degrees)

    cols = st.columns(2)
    with cols[0]:
        '#### Решения на различные моменты времени'
        fname = join(dir, f'a{a}_y{y}_p{p}_test{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    with cols[1]:
        '#### Закон сохранения массы'
        fname = join(dir, f'm_a{a}_y{y}_p{p}_test{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    