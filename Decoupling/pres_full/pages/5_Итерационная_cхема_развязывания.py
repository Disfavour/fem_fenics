import streamlit as st


r"""
# Итерационная cхема развязывания

Приближенное решение для $\varrho_{n+1},\ \bm u_{n+1}$ на $k$-й итерации обозначим $\varrho_{n+1}^{k},\ \bm u_{n+1}^{k}$,
причем используются начальные приближения на предыдущем слое по времени

$\begin{aligned}
\varrho_{n+1}^{0} = \varrho_{n}, \quad \bm u_{n+1}^{0} = \bm u_{n}
\end{aligned}$

Новое приближение на новом слое по времени получается на $K$-й итерации

$\begin{aligned}
\frac{\varrho_{n+1}^{k+1} - \varrho_{n}} {\tau} + \operatorname{div}(\varrho_{n+1}^{k+1} \bm u_{n+1}^k) = 0
\end{aligned}$

$\begin{aligned}
\frac{\varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1} - \varrho_{n} \bm u_{n}}{\tau } + \operatorname{div}(\bm u_{n+1}^{k} \otimes \varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1}) + \operatorname{grad} p(\varrho_{n+1}^{k+1}) = 0
\end{aligned}$

"""

with st.expander("Вариационная постановка"):
    r"""

Аналогично вариационной постановке линеаризованной схемы получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\int \limits_\Omega \frac{\varrho_{n+1}^{k+1} - \varrho_{n}} {\tau} \varrho_t \ dx
    - \int \limits_\Omega \varrho_{n+1}^{k+1} \bm u_{n+1}^k \cdot \grad \varrho_t \ dx
    = 0
    \\[0.5 cm]
    &\int \limits_\Omega \frac{\varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
    - \int \limits_\Omega \bm u_{n+1}^{k} \otimes \varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1} \vcentcolon \grad \bm u_t \ dx
    + \int \limits_\Omega \grad p(\varrho_{n+1}^{k+1}) \cdot \bm u_t \ dx
    = 0
\end{aligned} 
$$
"""
