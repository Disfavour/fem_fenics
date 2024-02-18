import streamlit as st

r'''
## Вычислительный алгоритм

### Аппроксимация по времени

Равномерная сетка по времени с шагом $\tau$
$$
\varphi^n = \varphi(t_n);\quad t_n = n \tau;\quad n = 0, 1, \dots, N;\quad N \tau = T
$$

Аппроксимация разностной производной
$$
{\left( \frac {\partial \varphi} {\partial t} \right)}^{n+1} \approx \frac {\varphi^{n+1} - \varphi^{n}} {\tau}
$$

Взвешенная аппроксимация
$$
\varphi^{n + \sigma} = \sigma \varphi^{n+1} + (1 - \sigma) \varphi^n, \quad \sigma = \operatorname{const}
$$

Схема с весами
$$
\begin{aligned}
	& \frac {h^{n+1} - h^n} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} u^{n+\sigma}) = 0
	\\[0.5 cm]
	& \frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} (u^{n+\sigma})^2) + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$

### Аппроксимация по пространству

Аппроксимация по пространству обеспечивается стандартными линейными конечными элементами
'''
