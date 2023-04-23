import streamlit as st


r"""
# Модельная задача

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

Параметры можно выбрать следующим образом:

$$
\begin{aligned}
&f = -6
\\[0.5 cm]
&\kappa_0 = 1
\\[0.5 cm]
&\kappa_1 = 0.01
\\[0.5 cm]
&u_D = 0
\end{aligned}
$$
"""
