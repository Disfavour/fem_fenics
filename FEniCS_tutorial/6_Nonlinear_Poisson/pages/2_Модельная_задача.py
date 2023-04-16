import streamlit as st


r"""
# Модельная задача



$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (q(u) \grad u) = f     \quad &x& \in \Omega,
\\[0.5 cm]
&u = u_D                         \quad &x& \in \partial \Omega,
\\[0.5 cm]
&\Omega = [0, 1] \times [0, 1].
\end{aligned}
$$

Пусть $u_e = 1 + x_1 + 2 x_2$ является решением, тогда остальные параметры можно выбрать следующим образом:

$$
\begin{aligned}
&q(u) = 1 + u^2
\\[0.5 cm]
&f = -10 - 10 x_1 - 20 x_2
\\[0.5 cm]
&u_D = u_e
\end{aligned}
$$
"""
