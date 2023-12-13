import streamlit as st

r'''
## Тестовая задача

Рассматривается двумерная задачи с возмущением плотности первоначально покоящейся жидкости

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	&\frac{\partial \varrho}{\partial t} + \div(\varrho \bm u) = 0 ,
	\quad &x& \in \Omega ,
	\quad &0& < t \leq T
	\\[0.5 cm]
	&\frac{\partial }{\partial t} (\varrho \bm u) + \div(\varrho \bm u \otimes \bm u) + \grad p = 0 , 
	\quad &x& \in \Omega ,
	\quad &0& < t \leq T
\end{aligned}
$$

Область
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	&\Omega = \{x \ | \ x = (x_1, x_2), \quad -5 < x_1 < 5, \quad -5 < x_2 < 5\}
\end{aligned}
$$

Начальные условия для скорости и плотности

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	&\bm u (x, 0) = \bm u^0(x) = 0
	\\[0.5 cm]
	&\varrho(x, 0) = \varrho^0(x) = 1 + \alpha \operatorname{exp}(-\beta |x|^2)
\end{aligned}
$$

Параметры

$$
\alpha = 2
\\[0.5 cm]
\beta = 20
$$



Давление определено следующим образом

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    & p(\varrho)  = a\varrho^\gamma
\end{aligned}
$$

Параметры

$$
a = 1
\\[0.5 cm]
\gamma = 2
$$
'''
