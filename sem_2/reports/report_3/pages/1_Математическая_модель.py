import streamlit as st

r'''
## Математическая модель

Уравнение неразрывности
$$
\frac {\partial h} {\partial t} + \operatorname{div} (h \bm u) = 0
$$

Уравнение движения в консервативном виде
$$
\frac {\partial} {\partial t} (h \bm u) + \operatorname{div} (h \bm u \otimes \bm u) + g h \operatorname{grad} h = 0
$$

### Одномерный случай

Уравнения мелкой воды
$$
\begin{aligned}
& \frac {\partial h} {\partial t} + \frac {\partial} {\partial x} (h u) = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (h u) + \frac {\partial} {\partial x} (h u^2) + g h \frac {\partial h} {\partial x} = 0
\end{aligned}
$$

Начальные условия
$$
\begin{aligned}
& h (x, 0) = h_0 (x)
\\[0.5cm]
& u (x, 0) = u_0 (x)
\end{aligned}
$$

Полная механическая энергия
$$
E =  \int_{\Omega} \frac{1}{2} \left( h u^2 + g h^2 \right) \ dx
$$
'''
