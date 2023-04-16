import streamlit as st


r"""
# Уравнение теплопроводности

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {\partial u} {\partial t} = \div \grad u + f         \quad &x& \in \Omega \times (0, T],
\\[0.5 cm]
&u = u_D                         \quad &x& \in \partial \Omega (0, T],
\\[0.5 cm]
&u = u_0 \quad &t& = 0.
\end{aligned}
$$

Здесь $u$ изменяется в зависимости от пространства и времени ( $u = u(x_1, x_2, \dots, t)$). Исходная функция $f$ и
граничные значения $u_D$ также могут изменяться в зависимости от пространства и времени. Начальное условие $u_0$
является функцией только пространства.

"""
