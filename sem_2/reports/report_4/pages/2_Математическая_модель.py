import streamlit as st

r'''
## Математическая модель

Координатная форма двумерных уравнений мелкой воды:
$$
	\begin{cases}
		\displaystyle\frac {\partial h} {\partial t}
		+ \frac {\partial} {\partial x_1} (h u)
		+ \frac {\partial} {\partial x_2} (h v) = 0
		\\[0.5 cm]
		\displaystyle\frac {\partial} {\partial t} (h u)
		+ \frac {\partial} {\partial x_1} (h u^2)
		+ \frac {\partial} {\partial x_2} (h u v)
		+ g h \frac {\partial h} {\partial x_1} = 0
		\\[0.5 cm]
		\displaystyle\frac {\partial} {\partial t} (h v)
		+ \frac {\partial} {\partial x_1} (h u v)
		+ \frac {\partial} {\partial x_2} (h v^2)
		+ g h \frac {\partial h} {\partial x_2} = 0
	\end{cases},
$$
где $g$ - ускорение свободного падения, $h(\bm x, t) > 0$ - высота слоя жидкости, $u(\bm x, t)$ - компонента горизонтальной скорости $\bm u (\bm x, t)$ по оси $x_1$, $v (\bm x, t)$ - компонента скорости по оси $x_2$.

Консервативная форма уравнений мелкой воды включает уравнение неразрывности:
$$
\begin{equation}
	\frac {\partial h} {\partial t} + \operatorname{div} (h \bm u) = 0
\end{equation}
$$
и уравнение движения:
$$
\begin{equation}
	\frac {\partial} {\partial t} (h \bm u) + \operatorname{div} (h \bm u \otimes \bm u) + g h \operatorname{grad} h = 0,
\end{equation}
$$
где символом $\otimes$ обозначено внешнее произведение (тензорное произведение двух векторов).

Боковая граница считается твердой. В силу этого имеем граничное условие непротекания:
$$
\begin{equation}
(\bm u \cdot \bm n) = 0, \quad \bm x \in \partial \Omega.
\end{equation}
$$

Задаются также начальные условия для высоты слоя жидкости и для скорости:
$$
\begin{equation}
h(\bm x, 0) = h_0 (\bm x), \quad \bm u (\bm x, 0) = \bm u_0 (\bm x), \quad \bm x \in \Omega.
\end{equation}
$$

Начально-краевая задача (1) - (4) описывает двумерные нестационарные течения несжимаемой жидкости
со свободной границей в приближении мелкой воды.

Непосредственным интегрированием уравнения неразрывности (1) по области $\Omega$ с учетом граничного
условия (3) получаем закон сохранения массы:
$$
\begin{equation}
m(t) = m(0), \quad m(t) = \varrho \int_{\Omega} h(\bm x, t) \ d \bm x,
\end{equation}
$$
где $\varrho = \operatorname{const} > 0$ - плотность жидкости.

Уравнение (2) напрямую выражает закон сохранения импульса. Интегрируя это уравнение по $\Omega$,
получим
$$
\int_\Omega \frac {\partial} {\partial t} (h \bm u) \ d \bm x + \frac g 2 \int_{\partial \Omega} h^2 \bm n \ d \bm x = 0.
$$
Тем самым имеет место
$$
\begin{equation}
\bm I (t) = \bm I (0) - \frac {g \varrho} 2 \int_0^t \int_{\partial \Omega} h^2 \bm n \ d \bm x \ dt, \quad \bm I (t) = \varrho \int_\Omega h \bm u \ d \bm x.
\end{equation}
$$

Принимая во внимание уравнение (1) и равенство
$$
\frac 1 2 \frac {\partial} {\partial t} (h \bm u^2) = \bm u \frac {\partial} {\partial t} (h \bm u) - \frac 1 2 \bm u^2 \frac {\partial h} {\partial t}
$$
и домножая уравнение (2) на $\bm u$, получим
$$
\begin{equation}
\frac 1 2 \frac {\partial} {\partial t} (h \bm u^2) + \frac 1 2 \operatorname{div} (h \bm u^2 \bm u) + g h \bm u \operatorname{grad} h = 0.
\end{equation}
$$
Уравнение (1) домножим на $g h$:
$$
\begin{equation}
\frac g 2 \frac {\partial h^2} {\partial t} + g h \operatorname{div} (h \bm u) = 0.
\end{equation}
$$
Складывая (7), (8), приходим к равенству
$$
\frac 1 2 \frac {\partial} {\partial t} (h \bm u^2) + \frac g 2 \frac {\partial h^2} {\partial t}
+ \frac 1 2 \operatorname{div} (h \bm u^2 \bm u) + g \operatorname{div} (h^2 \bm u) = 0.
$$
Интегрирование по области $\Omega$ с учетом (3) дает
$$
\frac d {dt} \int_\Omega \left( \frac 1 2 h \bm u^2 + \frac g 2 h^2 \right) \ d \bm x = 0.
$$
Приходим к закону сохранения полной механической энергии
$$
\begin{equation}
E(t) = E(0), \quad E(t) = \varrho \int_\Omega \left( \frac 1 2 h \bm u^2 + \frac g 2 h^2 \right) \ d \bm x.
\end{equation}
$$

Равенства (5), (6) и (9) являются основными законами сохранения для задачи (1) - (4).
'''
