import streamlit as st
from PIL import Image
import os.path


image = Image.open(os.path.join(os.path.dirname(__file__), 'im.png'))

r"""
# Более общие модели

Учет силы Кориолиса и рельефа дна
$$
\begin{aligned}
& \frac {\partial h} {\partial t} + \operatorname{div} (a \bm u) = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (a u)
+ \frac {\partial} {\partial x_1} (a u^2)
+ \frac {\partial} {\partial x_2} (a u v)
- f a v + g a \frac {\partial h} {\partial x_1} + c_f u \sqrt{u^2 + v^2} = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (a v)
+ \frac {\partial} {\partial x_1} (a u v)
+ \frac {\partial} {\partial x_2} (a v^2)
+ f a u + g a \frac {\partial h} {\partial x_2} + c_f v \sqrt{u^2 + v^2} = 0
\end{aligned}
$$

$f$ - параметр силы Кориолиса

$c_f$ - коэффициент для учета ветра
"""

st.image(image)
