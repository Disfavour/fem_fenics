import streamlit as st


r"""
# Курсовая работа

Конечно-элементное моделирование 2-мерного течения баротропной жидкости.

## Постановка задачи

- $\Omega$ - 2-мерная ограниченная область
- $\varrho (\bm x, t) > 0$ - плотность
- $\bm u (\bm x, t)$ - скорость
- $p (\bm x, t)$ - давление

Уравнение неразрывности

$$
\frac {\partial \varrho} {\partial t} + \operatorname{div} (\varrho \bm u) = 0, \quad \bm x \in \Omega
$$

Уравнение движения (в консервативном форме)

$$
\frac {\partial} {\partial t} (\varrho \bm u) + \operatorname{div} (\varrho \bm u \otimes \bm u)
+ \operatorname{grad} p = 0, \quad \bm x \in \Omega
$$

Уравнение состояния (баротропная жидкость)

$$
p = p(\varrho), \quad \frac {d p} {d \varrho} > 0 \quad (p = a \varrho^\gamma)
$$
"""
