import streamlit as st


r"""
# Математическая модель

Уравнение неразрывности в ограниченной области $\Omega$ имеет вид

$\begin{aligned}
\frac {\partial \varrho} {\partial t} + \operatorname{div} (\varrho \bm u) = 0,\quad x \in \Omega,\quad 0 < t \leq T
\end{aligned}$

где $\varrho(x, t) > 0$ - плотность, а $\bm u (x, t)$ - скорость

Уравнение движения в консервативном виде

$\begin{aligned}
\frac{\partial} {\partial t} (\varrho \bm u) + \operatorname{div}(\varrho \bm u \otimes \bm u) + \operatorname{grad} p = 0,\quad x \in \Omega,\quad 0 < t \leq T
\end{aligned}$

где $p(x, t)$ - давление

Жидкость предполагается баротропной, так что предполагается известной зависимость давления от плотности:
$p = p(\varrho), \ {\displaystyle \frac{d p}{d \varrho} > 0}$

Границы считаются твердыми, поэтому имеем граничное условие непротекания

$\begin{aligned}
(\bm u \cdot \bm n) = 0, \quad x \in \partial \Omega
\end{aligned}$

Задаются также начальные условия для плотности и скорости

$\begin{aligned}
\varrho(x, 0) = \varrho^0(x), \quad \bm u (x, 0) = \bm u^0(x), \quad x \in \Omega
\end{aligned}$

Начально-краевая задача описывает нестационарные течения идеальной баротропной жидкости
"""
