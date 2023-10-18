import streamlit as st

r'''
## Законы сохранения массы

Интегрированием уравнения неразрывности по области $\Omega$ с учетом граничного условия можно 
получить закон сохранения массы
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx + \int \limits_\Omega \div(\varrho \bm u) \ dx = 0
\end{aligned}
$$

Применяем теоремы о дивергенции
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx
+ \int \limits_{\partial \Omega} \varrho \bm u \cdot n \ dx = 0
\end{aligned}
$$

Учтем граничное условие:
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx = 0
\\[0.5 cm]
&\frac {\partial} {\partial t} \int \limits_\Omega \varrho \ dx = 0
\end{aligned}
$$

Таким образом, масса не зависит от времени:
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \varrho(x, t) \ dx = \int \limits_\Omega \varrho(x, 0) \ dx = \int \limits_\Omega \varrho^0(x) \ dx
\\[0.5cm]
& m (t) = m(0)
\\[0.5cm]
& m(t) = \int_{\Omega}  \varrho(x, t) d x
\end{aligned}
$$
'''
