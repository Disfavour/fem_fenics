import streamlit as st


r"""
# Математическая модель

$$
\bm u = (u, v)
$$

$$
\begin{aligned}
& \frac {\partial \varrho} {\partial t} + \operatorname{div} (\varrho \bm u) = 0
\\[0.5 cm]
& \frac {\partial u} {\partial t} + u \frac {\partial u} {\partial x_1} + v \frac {\partial u} {\partial x_2} + \frac {1} {\varrho}
\frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \frac {\partial v} {\partial t} + u \frac {\partial v} {\partial x_1} + v \frac {\partial v} {\partial x_2} + \frac {1} {\varrho}
\frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$

## Граничные условия

Граничное условие непротекания (no-slip)
$$
\bm u = 0 \qquad x \in \partial \Omega
$$

Граничное условие свободного скольжения (free-slip)
$$
\frac {\partial \bm u} {\partial n} = 0 \qquad x \in \partial \Omega
$$

## Начальные условия
$$
\begin{aligned}
& \bm u(x, 0) = \bm u^0(x)
\\
& \varrho(x) = \varrho^0
\end{aligned}
$$

## Закон сохранения массы
$$
\frac {\partial} {\partial t} \int \limits_\Omega \varrho \ dx + \int \limits_{\partial \Omega} \varrho \bm u \cdot n \ dx = 0
$$
"""
