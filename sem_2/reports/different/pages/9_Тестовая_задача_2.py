import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

dir = join(dirname(dirname(__file__)), 'images', 'err_etc')
hl, hr = 2, 1

r'''
# Тестовая задача 2
$$
\begin{aligned}
	& \frac {\partial h} {\partial t} + \frac {\partial} {\partial x} (h u) = 0
    \quad &x& \in \Omega,
	\quad &0& < t \leq T
	\\[0.5 cm]
	& \frac {\partial} {\partial t} (h u) + \frac {\partial} {\partial x} (h u^2) + g h \frac {\partial h} {\partial x} = 0
    \quad &x& \in \Omega,
	\quad &0& < t \leq T
\end{aligned}
$$

Область
$$
\Omega = \{ x \ | -5 \le x \le 5 \}
$$

Граничные условия
$$
u(-5, t) = u(5, t) = 0
$$

Начальные условия
$$
\begin{aligned}
	& h (x, 0) =
    \begin{cases}
		2, \quad &-&5 &\le &x &\le &0	\\[0.3 cm]
		1,  \quad  &&0  &< &x &\le &5
	\end{cases}
	\\[0.5cm]
	& u (x, 0) = 0
\end{aligned}
$$

## Аппроксимация по времени
Пусть $\tau$ - шаг равномерной сетки по времени
$$
\begin{aligned}
	&t_n = n \tau, \ n = 0, 1, \dots, N
    \\[0.5cm]
    &N \tau = T
\end{aligned}
$$

Номер временного слоя обозначим верхним индексом
$$
\varphi^n = \varphi(t_n)
$$

Таким образом, производная по времени может быть аппроксимирована
$$
{\left( \frac {\partial y} {\partial t} \right)}^{n+1} \approx \frac {y^{n+1} - y^{n}} {\tau}
$$

## Схема с весами
$$
\varphi^{n + \sigma} = \sigma \varphi^{n+1} + (1 - \sigma) \varphi^n, \quad \sigma = \operatorname{const}
$$

В результате система принимает вид
$$
\begin{aligned}
	& \frac {h^{n+1} - h^n} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} u^{n+\sigma}) = 0
	\\[0.5 cm]
	& \frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} (u^{n+\sigma})^2) + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''

r'## Решения на различные моменты времени'
cols = st.columns(3)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=1)
with cols[1]:
    tau = st.select_slider(r'$\tau$', (0.01, 0.005, 0.0025, 0.00125), 0.005, key=2)
with cols[2]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=3)
cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_h_ms{ms}_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_u_ms{ms}_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

r'## Влияние плотности сетки ($M$) на полную механическую энергию и нормы ошибок решений'
cols = st.columns(2)
with cols[0]:
    tau = st.select_slider(r'$\tau$', (0.01, 0.005, 0.0025, 0.00125), 0.005, key=4)
with cols[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=5)
with st.columns([0.25, 0.5, 0.25])[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_ms_tau{tau}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

r'## Влияние шага по времени ($\tau$) на полную механическую энергию и нормы ошибок решений'
cols = st.columns(2)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=6)
with cols[1]:
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1.0, 1.25), 1.0, key=7)
with st.columns([0.25, 0.5, 0.25])[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
cols = st.columns(2)    
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_tau_ms{ms}_s{s}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')

r'## Влияние весового параметра ($\sigma$) на полную механическую энергию и нормы ошибок решений'
cols = st.columns(2)
with cols[0]:
    ms = st.select_slider(r'$M$', (100, 200, 400, 800), 200, key=8)
with cols[1]:
    tau = st.select_slider(r'$\tau$', (0.01, 0.005, 0.0025, 0.00125), 0.005, key=9)
with st.columns([0.25, 0.5, 0.25])[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_E_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
cols = st.columns(2)
with cols[0]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_h_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
with cols[1]:
    fname = join(dir, f'hl{hl}_hr{hr}_err_u_dif_s_ms{ms}_tau{tau}.png')
    if exists(fname):
        st.image(Image.open(fname))
    else:
        st.error('')
