import streamlit as st


st.set_page_config(layout="wide")


r"""
# Краевая задача для нелинейного уравнения Пуассона

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (q(u) \grad u) = f     \quad &\bm{x}& \in \Omega
\\[0.5 cm]
&u = u_D                         \quad &\bm{x}& \in \partial \Omega
\end{aligned}
$$

"""
