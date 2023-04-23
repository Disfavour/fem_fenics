import streamlit as st


r"""
# Модельная задача

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {\partial u} {\partial t} = \div \grad u + f,         \quad &x& \in \Omega, \ &t \in (0, T]
\\[0.5 cm]
&u = u_D,                         \quad &x& \in \partial \Omega, \ &t \in (0, T]
\\[0.5 cm]
&u = u_0, \quad &t& = 0
\\[0.5 cm]
&\Omega = [-2, \ 2] \times [-2, \ 2]
\end{aligned}
$$

Параметры можно выбрать следующим образом:

$$
\begin{aligned}
&u_0(x_1, \ x_2) = e^{-a x_1^2 - a x_2^2}
\\[0.5 cm]
&u_D = 0
\\[0.5 cm]
&f = 0
\\[0.5 cm]
&a = 5
\end{aligned}
$$
"""
