import streamlit as st
import imageio


r"""
# Двумерный случай
"""

im = imageio.imread_v2("FEniCS_tutorial/4_function_space/d2.png")
st.image(im)

r"""
Двумерный случай отличается тем, что область разбивается на конечное число треугольников, которые являются конечными
элементами.

Аппроксимация функции на конечном элементе по 3-м узлам (вершины треугольника):

$$
\varphi(x) = \alpha_1 + \alpha_2 x_1 + \alpha_3 x_2
$$

## Лагранжевы элементы более высокой степени

Для элементов более высокой степени задействуются точки внутри и на ребрах треугольника:
"""

im = imageio.imread_v2("FEniCS_tutorial/4_function_space/d2m.png")
st.image(im)
