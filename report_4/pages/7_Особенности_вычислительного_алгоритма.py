import streamlit as st


r"""
# Особенности вычислительного алгоритма

## Новые переменные

$$
\begin{aligned}
& s = \varrho^\frac {1} {2}
\\[0.5 cm]
& \bm w = \varrho^\frac {1} {2} \bm u
\end{aligned}
$$

## Разностная схема Кранка—Николсона (симметричная схема)
$$
\begin{aligned}
& \frac {\partial u} {\partial t} = \frac {\partial^2 y} {\partial x^2}
\\[0.5 cm]
& \frac {\partial^2 y} {\partial x^2} = \frac {1} {2} \frac {\partial^2 u} {\partial x^2}
+ \frac {1} {2} \frac {\partial^2 u} {\partial x^2}
\\[0.5 cm]
& \frac {u^{n+1} - u^n} {\tau} = \frac {1} {2} \frac {\partial^2 u^{n+1}} {\partial x^2}
+ \frac {1} {2} \frac {\partial^2 u^{n}} {\partial x^2}
= \frac {\partial^2 u^{n+\frac {1} {2}}} {\partial x^2}
\end{aligned}
$$

## Законы сохранения на дискретном уровне

- Закон сохранения массы
- Закон сохранения импульса
- Закон сохранения энергии

"""
