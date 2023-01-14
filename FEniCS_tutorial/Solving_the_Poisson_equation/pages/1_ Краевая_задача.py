import streamlit as st


r"""
# Краевая задача Дирихле для уравнения Пуассона

Уравнение Пуассона:

$$
\begin{equation}
- \operatorname{div} \operatorname{grad} u = f(\bm x), \quad \bm x \in \Omega.
\end{equation}
$$

Граничные условия Дирихле:

$$
\begin{equation}
u(\bm x) = g(\bm x), \quad \bm x \in \partial \Omega.
\end{equation}
$$

Уравнение Пуассона в двумерном случае:

$$
\begin{equation}
- \frac {\partial^2 u} {\partial x_1^2} - \frac {\partial^2 u} {\partial x_2^2} = f(x_1, x_2).
\end{equation}
$$
"""
