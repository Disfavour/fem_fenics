import streamlit as st


r"""
# Линеаризованная схема развязывания (decoupling)

Из линейного уравнения может быть найдена плотность на новом слое по времени, затем из линейного векторного уравнения
может быть получена скорость

$\begin{aligned}
\frac {\varrho_{n+1} - \varrho_{n}} {\tau} + \operatorname{div}(\varrho_{n+1} \bm u_{n}) = 0
\end{aligned}$

$\begin{aligned}
\frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} + \operatorname{div}(\bm u_{n} \otimes \varrho_{n+1} \bm u_{n+1}) + \operatorname{grad} p(\varrho_{n+1}) = 0
\end{aligned}$
"""

with st.expander("Вариационная постановка"):
    r"""

## Вариационная постановка $(\varrho_{n+1})$

Умножаем на тестовую функцию и интегрируем по области

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
    + \int \limits_\Omega \div(\varrho_{n+1} \bm u_{n}) \varrho_t \ dx
    = 0
\end{aligned}
$$

Применяем формулу

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	& \div(\varphi \bm g) = \bm g \cdot \grad \varphi + \varphi \div \bm g
\end{aligned}
$$

Получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
    + \int \limits_\Omega \div(\varrho_{n+1} \varrho_t \bm u_{n}) \ dx
    - \int \limits_\Omega \varrho_{n+1} \bm u_{n} \cdot \grad \varrho_t \ dx
    = 0
\end{aligned}
$$

Применяем теорему о дивергенции

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
    + \int \limits_\Omega \varrho_{n+1} \varrho_t \bm u_{n} \cdot \bm n \ ds
    - \int \limits_\Omega \varrho_{n+1} \bm u_{n} \cdot \grad \varrho_t \ dx
    = 0
\end{aligned}
$$

С учетом граничных условий имеем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
    - \int \limits_\Omega \varrho_{n+1} \bm u_{n} \cdot \grad \varrho_t \ dx
    = 0
\end{aligned}
$$

## Вариационная постановка $(u_{n+1})$

Аналогично получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    \int \limits_\Omega \frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
    - \int \limits_\Omega \bm u_{n} \otimes \varrho_{n+1} \bm u_{n+1} \vcentcolon \grad \bm u_t \ dx
    + \int \limits_\Omega \grad p(\varrho_{n+1}) \cdot \bm u_t \ dx
    = 0
\end{aligned}
$$
"""
