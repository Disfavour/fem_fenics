import streamlit as st

r'''
## Схема с весами
$$
\varphi = \sigma \varphi^{n+1} + (1 - \sigma) \varphi^n, \quad \sigma = \operatorname{const}
$$

Пусть $\tau$ - шаг равномерной сетки во времени, такой, что
$$
\varphi^n = \varphi(t_n);\quad t_n = n \tau;\quad n = 0, 1, . . ., N;\quad N \tau = T
$$

Тогда
$$
\frac {\partial \varphi} {\partial t} = \frac {\varphi^{n+1} - \varphi^n} {\tau}
$$

Введем обозначение
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\def \um {(\sigma \bm u^{n+1} + (1 - \sigma) \bm u^{n})}
\begin{aligned}
&\bm u_w = \um
\\[0.5 cm]
&\varrho_w = (\sigma \varrho^{n+1} + (1 - \sigma) \varrho^{n})
\end{aligned}
$$

Решение на новом слое определяется из системы уравнений
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\def \um {(\sigma \bm u^{n+1} + (1 - \sigma) \bm u^{n})}
\begin{aligned}
& \frac {\varrho^{n+1} - \varrho^{n}} {\tau} + \div (\varrho_w \bm u_w) = 0
\\[0.5 cm]
& \frac {(\varrho \bm u)^{n+1} - (\varrho \bm u)^{n}} {\tau} + \div (\varrho_w \bm u_w \otimes \bm u_w)
+ \grad p = 0
\end{aligned}
$$
'''
