import streamlit as st


st.set_page_config(layout="wide")


r"""
# Краевая задача для нелинейного эллиптического уравнения 2-го порядка

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (q(u) \grad u) = f     \quad &x& \in \Omega,
\\[0.5 cm]
&u = u_D                         \quad &x& \in \partial \Omega.
\end{aligned}
$$
"""
