import streamlit as st


r'''
## Начально-краевая задача

Одномерные уравнения мелкой воды
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
\Omega = \{ x \ | -5 \le x \le 5 \}
$$

Граничные условия
$$
u(-5, t) = u(5, t) = 0
$$

Начальные условия
$$
\begin{aligned}
	& h (x, 0) =
    \begin{cases}
		h_l, \quad &-&5 &\le &x &\le &0	\\[0.3 cm]
		1,  \quad  &&0  &< &x &\le &5
	\end{cases}
	\\[0.5cm]
	& u (x, 0) = 0
\end{aligned}
$$
Рассматриваются задачи с разным параметром $h_l$
'''
