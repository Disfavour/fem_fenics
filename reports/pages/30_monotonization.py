import streamlit as st


r'''
## Монотонизация
$$
\frac {\partial h} {\partial t} + \frac \partial {\partial x} (h u) = 0
$$

#### Монотонизация 1
$$
\begin{aligned}
    &\frac {\partial h} {\partial t}
    - a \tau^2 \operatorname{div} (b \operatorname{grad} \frac {\partial h} {\partial t})
    + \frac \partial {\partial x} (h u) = 0
\end{aligned}
$$

Преобразования для вязкости в вариационной постановке
$$
\begin{aligned}
    &-\int \limits_\Omega a \tau^2 \operatorname{div} (b \operatorname{grad} \frac {\partial h} {\partial t}) v \ dx
    \\[0.5 cm]
    &a \tau^2 \left(
        -\int \limits_\Omega \operatorname{div} (b v \operatorname{grad} \frac {\partial h} {\partial t}) \ dx
        +\int \limits_\Omega b \operatorname{grad} \frac {\partial h} {\partial t} \operatorname{grad} v \ dx
    \right)
    \\[0.5 cm]
    &a \tau^2 \int \limits_\Omega b \operatorname{grad} \frac {\partial h} {\partial t} \operatorname{grad} v \ dx
\end{aligned}
$$

#### Монотонизация 2
$$
\begin{aligned}
    &\frac {\partial h} {\partial t}
    - a \tau^2 \operatorname{div} (b \operatorname{grad} h)
    + \frac \partial {\partial x} (h u) = 0
\end{aligned}
$$

Преобразования для вязкости в вариационной постановке
$$
\begin{aligned}
    &-\int \limits_\Omega a \tau^2 \operatorname{div} (b \operatorname{grad} h) v \ dx
    \\[0.5 cm]
    &a \tau^2 \left(
        -\int \limits_\Omega \operatorname{div} (b v \operatorname{grad} h) \ dx
        +\int \limits_\Omega b \operatorname{grad} h \operatorname{grad} v \ dx
    \right)
    \\[0.5 cm]
    &a \tau^2 \int \limits_\Omega b \operatorname{grad} h \operatorname{grad} v \ dx
\end{aligned}
$$

#### Замечания

Следующие выражения при вычислении в fenics на каждой итерации = 0
$$
\begin{aligned}
    &- a \tau^2 \operatorname{div} (b \operatorname{grad} \frac {\partial h} {\partial t})
    \\[0.5 cm]
    &- a \tau^2 \operatorname{div} (b \operatorname{grad} h)
\end{aligned}
$$

После обнуления слагаемого становятся некоторыми значениями
'''
