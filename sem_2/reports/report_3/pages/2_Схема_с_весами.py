import streamlit as st

r'''
## Схема с весами
$$
\varphi^{n + \sigma} = \sigma \varphi^{n+1} + (1 - \sigma) \varphi^n
$$

Пусть $\tau$ - шаг равномерной сетки во времени, такой, что
$$
\varphi^n = \varphi(t_n);\quad t_n = n \tau;\quad n = 0, 1, . . ., N;\quad N \tau = T
$$

Тогда
$$
\frac {\partial \varphi} {\partial t} = \frac {\varphi^{n+1} - \varphi^n} {\tau}
$$

Решение на новом временном слое определяется из системы уравнений
$$
\def \s {{n+\sigma}}
\begin{aligned}
& \frac {h^{n+1} - h^n} {\tau} + \frac {\partial} {\partial x} (h^\s u^\s) = 0
\\[0.5 cm]
& \frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} + \frac {\partial} {\partial x} ( h^\s \displaystyle (u^\s)^2) + g h^\s \frac {\partial h^\s} {\partial x} = 0
\end{aligned}
$$
'''
