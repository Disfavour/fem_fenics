import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


name_ = 'weight'
mesh_size = 200
tau = 0.0025
tests = list(range(1, 4))
degrees = list(range(1, 4))
kappas = [0] + [0.01 * 2**i for i in range(8)]

dir = join(dirname(dirname(dirname(dirname(__file__)))), 'images', 'monotonization', 'linear', name_)

'## Монотонизирующая схема с весами'

name = st.sidebar.radio("Содержание", (
    "Вычислительная схема",
    'Тестовые задачи',
    "Численные эксперименты",
    ))

if name == "Вычислительная схема":
    r'''
Уравнение переноса
$$
\frac {\partial u} {\partial t} + \frac {\partial u} {\partial x} = 0
$$

Симмертичная двухслойная схема
$$
\frac {u^{n+1} - u^n} \tau + \frac \partial {\partial x} \frac {u^{n+1} + u^n} 2 = 0
$$

Операторный вид
$$
\frac {u^{n+1} - u^n} \tau + A \frac {u^{n+1} + u^n} 2 = 0
$$

Добавленное слагаемое
$$
\left( 1 + \left( \sigma - \frac 1 2 \right) \tau A \right) \frac {u^{n+1} - u^n} \tau + A \frac {u^{n+1} + u^n} 2 = 0
$$

При $\sigma = \displaystyle \frac 1 2$ полностью соответствует неявной схеме

Повышение веса повышает монотонность (аналогично весу двухслойной схемы с весами)

Замена $\sigma - \displaystyle\frac 1 2$ на $\kappa$
$$
\left( 1 + \kappa \tau A \right) \frac {u^{n+1} - u^n} \tau + A \frac {u^{n+1} + u^n} 2 = 0
$$
'''

elif name == 'Тестовые задачи':
    '### Тестовые задачи'
    r'''
$$
\frac {\partial u} {\partial t} + \frac {\partial u} {\partial x} = 0
$$

Область
$$
\Omega = \{ x \ | 0 \le x \le 1 \}
$$

Граничные условия
$$
u(0, t) = 0
$$

Начальные условия
$$
\begin{aligned}
	& u (x, 0) = u_i(x)
    \\[0.5 cm]
    &u_1(x) = 
    \begin{cases}
		1, \quad x \in [0.15, 0.25]	\\[0.3 cm]
		0, \quad x \notin [0.15, 0.25]
	\end{cases}
    \\[0.5 cm]
    &u_2(x) = 
    \begin{cases}
		1, \quad \sin(10 \pi (x - 0.15)) \in [0.15, 0.25]	\\[0.3 cm]
		0, \quad x \notin [0.15, 0.25]
	\end{cases}
    \\[0.5 cm]
    &u_3(x) = e^{-2000(x - 0.2)^2}
\end{aligned}
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
    