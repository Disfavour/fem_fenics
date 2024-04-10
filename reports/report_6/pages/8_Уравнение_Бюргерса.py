import streamlit as st
from PIL import Image
from os.path import dirname, join, exists


mesh_size = 200
tau = 0.0025
sigma = 0.5
ntests = (1, 2)

dir_base = join(dirname(dirname(dirname(dirname(__file__)))), 'images', 'monotonization', 'nonlinear')

'## Уравнение Бюргерса'

name = st.sidebar.radio("Содержание", (
    "Вычислительная схема",
    'Тестовые задачи',
    "Численные эксперименты",
    ))

if name == "Вычислительная схема":
    r'''
$$
\frac {\partial u} {\partial t} + u \frac {\partial u} {\partial x} = 0
$$

Монотонизирующая схема с весами
$$
\left( 1 + \kappa \tau A \right) \frac {\partial u} {\partial t} + u \frac {\partial u} {\partial x} = 0
$$

Нелинейная вязкость
$$
\begin{aligned}
    &\frac {\partial u} {\partial t} + u \frac {\partial u} {\partial x}
    - a \tau^2 \frac \partial {\partial x} \left( b \frac {\partial u} {\partial x} \right)
    = 0
    \\[0.5 cm]
    & a = \operatorname{const} > 0
    \\[0.5 cm]
    & b = | \operatorname{grad} u |^{\gamma}
\end{aligned}
$$
'''

elif name == 'Тестовые задачи':
    '### Тестовые задачи'
    r'''
$$
\frac {\partial u} {\partial t} + u \frac {\partial u} {\partial x} = 0
$$

Область
$$
\Omega = \{ x \ | 0 \le x \le 1 \}
$$

Граничные условия
$$
u(0, t) = a_i
$$

Начальные условия
$$
\begin{aligned}
	& u (x, 0) =
    \begin{cases}
		a_i, \quad 0    \le x \le x^m_i	\\[0.3 cm]
		b_i, \quad x^m_i  < x \le 1
	\end{cases}
    \\[0.8 cm]
    &a_1, b_1, x^m_1 = 0.5, 1.5, 0.2
    \\[0.5 cm]
    &a_2, b_2, x^m_2 = 1.5, 0.5, 0.2
\end{aligned}
$$
    '''

elif name == "Численные эксперименты":
    "### Численные эксперименты"
    with st.columns([0.3, 0.4, 0.3])[1]:
        test = st.selectbox('Тестовая задача', ('1', '2'))
    cols = st.columns(2)

    with cols[0]:
        "#### Монотонизирующая схема"
        dir_weight = join(dir_base, 'weight')
        degrees_weight = list(range(1, 4))
        alfas_weight = [0]
        gammas_weight = [0]
        kappas_weight = [0] + [0.01 * 2**i for i in range(8)]

        a_weight = alfas_weight[0]
        y_weight = gammas_weight[0]

        cols1 = st.columns(2)
        with cols1[0]:
            k_weight = st.select_slider(r'$\kappa$', kappas_weight)
        with cols1[1]:
            p_weight = st.select_slider(r'$p$', degrees_weight)
        
        '#### Решения на различные моменты времени'
        fname = join(dir_weight, f'u_ms{mesh_size}_tau{tau}_d{p_weight}_alfa{a_weight}_gamma{y_weight}_kappa{k_weight}_ntest{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
        
        '#### $L_2$ норма ошибки'
        fname = join(dir_weight, f'err_ms{mesh_size}_tau{tau}_d{p_weight}_alfa{a_weight}_gamma{y_weight}_kappa{k_weight}_ntest{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
    
    with cols[1]:
        "#### Нелинейная вязкость"
        dir = join(dir_base, 'viscosity')
        degrees = list(range(1, 4))
        alfas = [0] + [2**i for i in range(8)]
        gammas = list(range(4))
        kappas = [0]

        k = kappas[0]

        cols1 = st.columns(3)
        with cols1[0]:
            a = st.select_slider(r'$a$', alfas)
        with cols1[1]:
            y = st.select_slider(r'$\gamma$', gammas)
        with cols1[2]:
            p = st.select_slider(r'$p$', degrees, key=123)
        
        '#### Решения на различные моменты времени'
        fname = join(dir, f'u_ms{mesh_size}_tau{tau}_d{p}_alfa{a}_gamma{y}_kappa{k}_ntest{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
        
        '#### $L_2$ норма ошибки'
        fname = join(dir, f'err_ms{mesh_size}_tau{tau}_d{p}_alfa{a}_gamma{y}_kappa{k}_ntest{test}.png')
        if exists(fname):
            st.image(Image.open(fname))
        else:
            st.error('')
