import streamlit as st

r'''
## Вариационная формулировка
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\def \um {(\sigma \bm u^{n+1} + (1 - \sigma) \bm u^{n})}
\begin{aligned}
& \frac {\varrho^{n+1} - \varrho^{n}} {\tau} + \div (\varrho_w \bm u_w) = 0
\\[0.5 cm]
& \frac {(\varrho \bm u)^{n+1} - (\varrho \bm u)^{n}} {\tau} + \div (\varrho_w \bm u_w \otimes \bm u_w)
+ \grad p(\varrho_w) = 0
\end{aligned}
$$

Умножаем каждое уравнение на тестовую функцию, интегрируем по области и суммируем уравнения
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\begin{aligned}
& \int \limits_\Omega \frac {\varrho^{n+1} - \varrho^{n}} {\tau} \varrho_t \ dx
+ \int \limits_\Omega \div(\varrho_w \bm u_w) \varrho_t \ dx
+ \int \limits_\Omega \frac {\varrho^{n+1} \bm u^{n+1} - \varrho^{n} \bm u^{n}} {\tau} \cdot \bm u_t \ dx
\\[0.5 cm]
&+ \int \limits_\Omega \div(\varrho_w \bm u_w \otimes \bm u_w) \cdot \bm u_t \ dx
+ \int \limits_\Omega \grad p(\varrho_w)\cdot \bm u_t \ dx
= 0
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
& \int \limits_\Omega \frac {\varrho^{n+1} - \varrho^{n}} {\tau} \varrho_t \ dx
+ \int \limits_\Omega \div(\varrho_w \varrho_t \bm u_w) \ dx
- \int \limits_\Omega \varrho_w \bm u_w \cdot \grad \varrho_t \ dx
\\[0.5 cm]
&+ \int \limits_\Omega \frac {\varrho^{n+1} \bm u^{n+1} - \varrho^{n} \bm u^{n}} {\tau} \cdot \bm u_t \ dx
+ \int \limits_\Omega \div(\varrho_w \bm u_w \otimes \bm u_w \cdot \bm u_t) \ dx
\\[0.5 cm]
&- \int \limits_\Omega \varrho_w \bm u_w \otimes \bm u_w \vcentcolon \grad \bm u_t  \ dx
+ \int \limits_\Omega \grad p(\varrho_w)\cdot \bm u_t \ dx
= 0
\end{aligned}
$$

Применяем теорему о дивергенции

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	\begin{split}
	    & \int \limits_\Omega \frac {\varrho^{n+1} - \varrho^{n}} {\tau} \varrho_t \ dx
	    + \int \limits_{\partial \Omega} \varrho_w \varrho_t \bm u_w \cdot \bm n \ ds
	    - \int \limits_\Omega \varrho_w \bm u_w \cdot \grad \varrho_t \ dx
	    \\[0.5 cm]
		&+ \int \limits_\Omega \frac {\varrho^{n+1} \bm u^{n+1} - \varrho^{n} \bm u^{n}} {\tau} \cdot \bm u_t \ dx
		+ \int \limits_{\partial \Omega} \varrho_w \bm u_w \otimes \bm u_w \cdot \bm u_t \cdot \bm n \ ds
		\\[0.5 cm]
		&- \int \limits_\Omega \varrho_w \bm u_w \otimes \bm u_w \vcentcolon \grad \bm u_t  \ dx
	    + \int \limits_\Omega \grad p(\varrho_w)\cdot \bm u_t \ dx
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
	    & \int \limits_\Omega \frac {\varrho^{n+1} - \varrho^{n}} {\tau} \varrho_t \ dx
	    - \int \limits_\Omega \varrho_w \bm u_w \cdot \grad \varrho_t \ dx
		+ \int \limits_\Omega \frac {\varrho^{n+1} \bm u^{n+1} - \varrho^{n} \bm u^{n}} {\tau} \cdot \bm u_t \ dx
		\\[0.5 cm]
		&- \int \limits_\Omega \varrho_w \bm u_w \otimes \bm u_w \vcentcolon \grad \bm u_t  \ dx
	    + \int \limits_\Omega \grad p(\varrho_w)\cdot \bm u_t \ dx
	    = 0
	\end{split}
\end{aligned}
$$
'''
