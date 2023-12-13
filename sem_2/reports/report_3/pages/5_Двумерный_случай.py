import streamlit as st
from os.path import join, dirname

images_dir = join(dirname(dirname(__file__)), 'images')

sw_2d = open(join(images_dir, 'sw_2d.mp4'), 'rb').read()
sw_2d_smooth = open(join(images_dir, 'sw_2d_smooth.mp4'), 'rb').read()

r'''
## Двумерный случай
$$
\begin{aligned}
& \frac {\partial h} {\partial t} + \operatorname{div} (h \bm u) = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (h \bm u) + \operatorname{div} (h \bm u \otimes \bm u) + g h \operatorname{grad} h = 0
\end{aligned}
$$

Область
$$
\Omega = \{ \bm x | \bm x = \{x_1, x_2\}, \quad -5 <= x_1 <= 5, \quad -1 <= x_2 <= 1 \}
$$

Граничные условия
$$
\bm u \cdot \bm n = 0, \quad \bm x \in \partial \Omega
$$
'''

cols = st.columns(2)
with cols[0]:
    r'''
    ### Тестовая задача 1

Начальные условия
$$
\begin{aligned}
	& h (\bm x, 0) =
    \begin{cases}
		10, \quad &-&5 &<= &x_0 &<= &0	\\[0.3 cm]
		1,  \quad  &&0  &< &x_0 &<= &5
	\end{cases}
	\\[0.5cm]
	& \bm u (\bm x, 0) = 0
\end{aligned}
$$'''
    st.video(sw_2d)
with cols[1]:
    r'''
    ### Тестовая задача 2
Начальные условия
$$
\begin{aligned}
	& h (\bm x, 0) =
    \begin{cases}
		3, 					\quad &-&5 &<= &x_0 &<= &0	\\[0.3 cm]
		1 + 2 e^{-20 x_0^2},	\quad  &&0  &< &x_0 &<= &5
	\end{cases}
	\\[0.5cm]
	& \bm u (\bm x, 0) = 0
\end{aligned}
$$'''
    st.video(sw_2d_smooth)
