import streamlit as st


r"""
# Модельная задача

Рассматривается двумерная задачи с возмущением плотности первоначально покоящейся жидкости

$\begin{aligned}
\frac{\partial \varrho}{\partial t} + \operatorname{div}(\varrho \bm u) = 0, \quad\quad \quad\quad\quad\quad\quad\quad\ \ &x& \in \Omega ,\quad &0& < t \leq T
\end{aligned}$

$\begin{aligned}
\frac{\partial }{\partial t} (\varrho \bm u) + \operatorname{div}(\varrho \bm u \otimes \bm u) + \operatorname{grad} p = 0 , \quad &x& \in \Omega ,\quad &0& < t \leq T
\end{aligned}$

Начальные условия для скорости и плотности

$\begin{aligned}
\bm u (x, 0)  = 0
\end{aligned}$

$\begin{aligned}
\varrho(x, 0) = 1 + 2 \operatorname{exp}(-20 |x|^2)
\end{aligned}$


Расчетная область

$\begin{aligned}
\Omega = \{x \ | \ x = (x_1, x_2), \quad -5 < x_1 < 5, \quad -5 < x_2 < 5\}
\end{aligned}$

Давление

$\begin{aligned}
p(\varrho)  = \varrho^{1.4}
\end{aligned}$
"""
