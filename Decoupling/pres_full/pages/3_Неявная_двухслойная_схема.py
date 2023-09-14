import streamlit as st


r"""
# Неявная двухслойная схема

Пусть $\tau$ - шаг равномерной сетки во времени, такой, что

$\phi_n = \phi(t_n);\quad t_n = n \tau;\quad n = 0, 1, . . ., N;\quad N \tau = T$

Решение на новом слое определяется из уравнений


$\begin{aligned}
\frac {\varrho_{n+1} - \varrho_{n}} {\tau} + \operatorname{div}(\varrho_{n+1} \bm u_{n+1}) = 0
\end{aligned}$

$\begin{aligned}
\frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}  }{\tau } + \operatorname{div}(\varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1}) + \operatorname{grad} p(\varrho_{n+1}) = 0
\end{aligned}$
"""

with st.expander("Вариационная постановка"):
    r"""
Умножаем каждое уравнение на тестовую функцию, интегрируем по области и суммируем уравнения

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	\begin{split}
	    & \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
	    + \int \limits_\Omega \div(\varrho_{n+1} \bm u_{n+1}) \varrho_t \ dx
		+ \int \limits_\Omega \frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
		\\[0.5 cm]
		&+ \int \limits_\Omega \div(\varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1}) \cdot \bm u_t \ dx
	    + \int \limits_\Omega \grad p(\varrho_{n+1})\cdot \bm u_t \ dx
	    = 0
	\end{split}
\end{aligned}
$$

Применяем формулы

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	& \div(\varphi \bm g) = \bm g \cdot \grad \varphi + \varphi \div \bm g
	\\[0.5 cm]
	& \div(\bm A \cdot \bm f) = \bm A \vcentcolon \grad \bm f + \div \bm A \cdot \bm f
\end{aligned}
$$

Получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	\begin{split}
	    & \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
	    + \int \limits_\Omega \div(\varrho_{n+1} \varrho_t \bm u_{n+1}) \ dx
	    - \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \cdot \grad \varrho_t \ dx
	    \\[0.5 cm]
		&+ \int \limits_\Omega \frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
		+ \int \limits_\Omega \div(\varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1} \cdot \bm u_t) \ dx
		\\[0.5 cm]
		&- \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1} \vcentcolon \grad \bm u_t  \ dx
	    + \int \limits_\Omega \grad p(\varrho_{n+1})\cdot \bm u_t \ dx
	    = 0
	\end{split}
\end{aligned}
$$

Применяем теорему о дивергенции

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	\begin{split}
	    & \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
	    + \int \limits_{\partial \Omega} \varrho_{n+1} \varrho_t \bm u_{n+1} \cdot \bm n \ ds
	    - \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \cdot \grad \varrho_t \ dx
	    \\[0.5 cm]
		&+ \int \limits_\Omega \frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
		+ \int \limits_{\partial \Omega} \varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1} \cdot \bm u_t \cdot \bm n \ ds
		\\[0.5 cm]
		&- \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1} \vcentcolon \grad \bm u_t  \ dx
	    + \int \limits_\Omega \grad p(\varrho_{n+1})\cdot \bm u_t \ dx
	    = 0
	\end{split}
\end{aligned}
$$

С учетом граничных условий имеем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	\begin{split}
	    & \int \limits_\Omega \frac {\varrho_{n+1} - \varrho_{n}} {\tau} \varrho_t \ dx
	    - \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \cdot \grad \varrho_t \ dx
		+ \int \limits_\Omega \frac {\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}} {\tau} \cdot \bm u_t \ dx
		\\[0.5 cm]
		&- \int \limits_\Omega \varrho_{n+1} \bm u_{n+1} \otimes \bm u_{n+1} \vcentcolon \grad \bm u_t  \ dx
	    + \int \limits_\Omega \grad p(\varrho_{n+1})\cdot \bm u_t \ dx
	    = 0
	\end{split}
\end{aligned}
$$
"""
