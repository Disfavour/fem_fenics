import streamlit as st

r'''
## Тестовые задачи

$$
\begin{aligned}
	& \frac {\partial h} {\partial t} + \frac {\partial} {\partial x} (h u) = 0
    \quad &x& \in \Omega,
	\quad &0& < t \leq T
	\\[0.5 cm]
	& \frac {\partial} {\partial t} (h u) + \frac {\partial} {\partial x} (h u^2) + g h \frac {\partial h} {\partial x} = 0
    \quad &x& \in \Omega,
	\quad &0& < t \leq T
\end{aligned}
$$

Область
$$
\Omega = \{ x \ | -5 < x < 5 \}
$$

Граничные условия
$$
u(-5, t) = u(5, t) = 0
$$

### Тестовая задача 1

Начальные условия
$$
\begin{aligned}
	& h (x, 0) =
    \begin{cases}
		10, \quad &-&5 &<= &x &<= &0	\\[0.3 cm]
		1,  \quad  &&0  &< &x &<= &5
	\end{cases}
	\\[0.5cm]
	& u (x, 0) = 0
\end{aligned}
$$

### Тестовая задача 2
Начальные условия
$$
\begin{aligned}
	& h (x, 0) =
    \begin{cases}
		3, 					\quad &-&5 &<= &x &<= &0	\\[0.3 cm]
		1 + 2 e^{-20 x^2},	\quad  &&0  &< &x &<= &5
	\end{cases}
	\\[0.5cm]
	& u (x, 0) = 0
\end{aligned}
$$
'''
