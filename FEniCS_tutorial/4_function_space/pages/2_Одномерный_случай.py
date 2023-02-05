import streamlit as st
import imageio


r"""
# Одномерный случай
"""

im = imageio.imread_v2("FEniCS_tutorial/4_function_space/d1.png")
st.image(im)

r"""
Область:

$$
\Omega = \{x \ | \ 0 \le x \le l\}
$$

Сетка:

$$
x_j = \alpha h
\\[0.3 cm]
j = 0, 1, \dots, n
\\[0.3 cm]
\alpha = 0, 1, \dots, n
\\[0.3 cm]
h = \frac {l} {n}
$$

Конечные элементы:

$$
\Omega_i = \{ x \  | \  x_{j-1} \le x \le x_j  \}
\\[0.3 cm]
i = 1, \dots, n
$$

Аппроксимация полиномом по 2-м узлам:

$$
\varphi = \alpha_1 + \alpha_2 x
$$

Коэффициенты $\alpha_1$ и $\alpha_2$ могут быть определены с помощью условий в узловых точках:

- $\varphi = \varPhi_k$ при $x = x_k$;
- $\varphi = \varPhi_{k+1}$ при $x = x_{k+1}$.

Эти условия приводят к системе двух уравнений:

$$
\begin{cases}
\varPhi_k = \alpha_1 + \alpha_2 x_k
\\
\varPhi_{k+1} = \alpha_1 + \alpha_2 x_{k+1}
\end{cases}
$$

Полиномиальная функция:

$$
\varphi = \left( \frac {x_{k+1} - x} {x_{k+1} - x_k} \right) \varPhi_k + \left( \frac {x - x_k} {x_{k+1} - x_k} \right)
\varPhi_{k+1}
$$

Значения функции равны единице в одном определённом узле и обращаются в ноль во всех других узлах.

$$
u(x) = \sum \limits_{i=1}^n c_i \varphi_i(x)
\\[0.3 cm]
x \in \Omega = \bigcup \limits_{i=1}^m \Omega_i
$$

Конечно-элементный базис - функции $\varphi_i(x), \  i = 1, 2, \dots, n$, которые в узлах:

$$
\varphi_i(x_j) = 
\begin{cases}
1, \  i = j
\\
0, \  i \ne j
\end{cases}
$$


## Лагранжевы элементы более высокой степени

Для элементов более высокой степени задействуются точки внутри отрезка:
"""

im = imageio.imread_v2("FEniCS_tutorial/4_function_space/d1m.png")
st.image(im)
