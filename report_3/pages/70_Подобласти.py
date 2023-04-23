import streamlit as st
from PIL import Image
import os.path


image = Image.open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "images", 'subdomains.png'))


r"""
# Подобласти

Решение краевых задач в областях, состоящих из различных материалов, в FEniCS решаются путем определения
подобластей внутри области.

Рассмотрим следующую краевую задачу для нелинейного эллиптического уравнения:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (\kappa(x) \grad u) = f,        \quad &x& \in \Omega
\\[0.5 cm]
&u = u_D,                               \quad &x& \in \partial \Omega
\\[0.5 cm]
&\Omega_0 = [0,\ 1] \times [ 0,\ 0.5 ]
\\[0.5 cm]
&\Omega_1 = [0,\ 1] \times [ 0.5,\ 1 ]
\\[0.5 cm]
&\kappa(x) =
\begin{cases}
   \kappa_0, &x \in \Omega_0 \\
   \kappa_1, &x \in \Omega_1
\end{cases}
\end{aligned}
$$
"""

st.image(image, width=300)
