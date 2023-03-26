import streamlit as st


r"""
# Постановка задачи

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (q(u) \grad u) = f     \quad &\bm{x}& \in \Omega
\\[0.5 cm]
&u = u_D                         \quad &\bm{x}& \in \partial \Omega
\\[1.0 cm]
&\Omega = [0, 1] \times [0, 1]
\\[0.5 cm]
&q(u) = 1 + u^2
\\[0.5 cm]
&u_D = 1 + x_1 + 2 x_2
\\[0.5 cm]
&f = -10 - 10 x_1 - 20 x_2
\\[0.5 cm]
&u_e = 1 + x_1 + 2 x_2
\end{aligned}
$$
"""
