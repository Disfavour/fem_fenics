import streamlit as st


r"""
# Итерационная cхема развязывания

На основе линеаризованной схемы можно строить итерационный алгоритм для вычислительной реализации чисто неявной схемы

Приближенное решение для $\varrho_{n+1},\ \bm u_{n+1}$ на $k$-й итерации обозначим $\varrho_{n+1}^{k},\ \bm u_{n+1}^{k}$,
причем используется начальное приближение на предыдущем слое по времени

$$
\varrho_{n+1}^{0} = \varrho_{n}, \quad \bm u_{n+1}^{0} = \bm u_{n}
$$

Новое приближение на новом слое по времени находится, когда делается $K$ итераций

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \frac{\varrho_{n+1}^{k+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1}^{k+1} \bm u_{n+1}^k)
    = 0
\end{aligned} 
$$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \frac{\varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1} - \varrho_{n} \bm u_{n}}{\tau }
    + \div(\bm u_{n+1}^{k} \otimes \varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1})
    + \grad p(\varrho_{n+1}^{k+1})
    = 0
\end{aligned} 
$$

## Вариационная постановка

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
