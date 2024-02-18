import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

images = join(dirname(dirname(__file__)), 'images')

r'''
## Время и сигма

$$
\varphi^{n+\sigma} = \sigma \varphi^{n+1} + (1 - \sigma) \varphi^n
\\[0.5 cm]
\varphi^{n+1} = \frac {\varphi^{n+\sigma} - (1 - \sigma) \varphi^n} \sigma
\\[0.5 cm]
\frac {\partial \varphi} {\partial t} = \frac {\varphi^{n+1} - \varphi^n} \tau
$$

### Тестовая задача

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
\Omega = \{ x \ | -5 < x < 5 \}
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
		3, \quad &-&5 &\le &x &\le &0	\\[0.3 cm]
		1,  \quad  &&0  &< &x &\le &5
	\end{cases}
	\\[0.5cm]
	& u (x, 0) = 0
\end{aligned}
$$
'''

cols = st.columns(2)

with cols[0]:
    r'''### 1
$$
\begin{aligned}
	& \frac {h^{n+1} - h^n} {\tau} + \frac {\partial} {\partial x} (h u)^{n+\sigma} = 0
	\\[0.5 cm]
	& \frac {(h u)^{n+1} - (h u)^n} {\tau} + \frac {\partial} {\partial x} (h u^2)^{n+\sigma} + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''

with cols[1]:
    r'''### 2
$$
\begin{aligned}
	& \frac {h^{n+1} - h^n} {\tau} + \frac {\partial} {\partial x} (h u)^{n+\sigma} = 0
	\\[0.5 cm]
	& h^{n+\sigma} \frac {u^{n+1} - u^n} \tau + u^{n+\sigma} \frac {h^{n+1} - h^n} \tau + \frac {\partial} {\partial x} (h u^2)^{n+\sigma} + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''

cols = st.columns(2)

with cols[0]:
    r'''### 3
$$
\begin{aligned}
	& \frac {h^{n+\sigma} - h^n} {\tau \sigma} + \frac {\partial} {\partial x} (h u)^{n+\sigma} = 0
	\\[0.5 cm]
	& \frac {\frac {h^{n+\sigma} - (1 - \sigma) h^n} \sigma \frac {u^{n+\sigma} - (1 - \sigma) u^n} \sigma - (h u)^n} {\tau} + \frac {\partial} {\partial x} (h u^2)^{n+\sigma} + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''

with cols[1]:
    r'''### 4
$$
\begin{aligned}
	& \frac {h^{n+\sigma} - h^n} {\tau \sigma} + \frac {\partial} {\partial x} (h u)^{n+\sigma} = 0
	\\[0.5 cm]
	& h^{n+\sigma} \frac {\frac {u^{n+\sigma} - (1 - \sigma) u^n} \sigma - u^n} \tau + u^{n+\sigma} \frac {\frac {h^{n+\sigma} - (1 - \sigma) h^n} \sigma - h^n} \tau + \frac {\partial} {\partial x} (h u^2)^{n+\sigma} + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''

with st.columns(3)[1]:
    m = st.select_slider(r'$m$', (100, 200, 400), 200)
    tau = st.select_slider(r'$\tau$', (0.0025, 0.005, 0.01), 0.005)
    s = st.select_slider(r'$\sigma$', (0.5, 0.75, 1, 1.5), 1)

# ms100_tau0.01_s0.75_c0.png
cols = st.columns(2)
with cols[0]:
    fname = join(images, f'240107_ms{m}_tau{tau}_s{s}_E.png')
    if exists(fname):
            st.image(Image.open(fname))
with cols[1]:
    fname = join(images, f'240107_ms{m}_tau{tau}_s{s}_dif_E.png')
    if exists(fname):
            st.image(Image.open(fname))