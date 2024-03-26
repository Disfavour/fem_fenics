import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

images = join(dirname(dirname(__file__)), 'images')

r'''
## Уравнения с давлением

Базовая версия
$$
\begin{aligned}
& \frac {\partial h} {\partial t} + \operatorname{div} (h \bm u) = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (h \bm u) + \operatorname{div} (h \bm u \otimes \bm u) + g h \operatorname{grad} h = 0
\end{aligned}
$$

Давление
$$
p = a \varrho^\gamma
$$

Потенциал давления
$$
\varrho \frac {\partial \varPi} {\partial \varrho} - \varPi (\varrho) = p (\varrho)
$$

Получаем 
$$
\varPi (\varrho) = a \frac {\varrho^\gamma} {\gamma - 1}
$$

Уравнения неразрывности относительно давления:
$$
\frac {\partial \varPi} {\partial t} + \operatorname{div} (\varPi \bm u) + p(\varrho) \operatorname{div} \bm u = 0
$$

Рассматривается в приближении мелкой воды, поэтому
$$
a = \frac g 2
\\[0.5 cm]
\gamma = 2
\\[0.5 cm]
\varrho = h
$$

Выражаем $h$ через давление:
$$
h = \left( \frac p a \right)^{\displaystyle\frac 1 \gamma}
$$

### Версия 2
$$
\def \h {\left( \displaystyle\frac p a \right)^{\displaystyle\frac 1 \gamma}}
\begin{aligned}
& \frac {\partial p} {\partial t} + \operatorname{div} (p \bm u) + p \operatorname{div} \bm u = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (\h \bm u) + \operatorname{div} (\h \bm u \otimes \bm u) + \operatorname{grad} p = 0
\end{aligned}
$$

### Версия 3
$$
\def \h {\left( \displaystyle\frac p a \right)^{\displaystyle\frac 1 \gamma}}
\begin{aligned}
& \frac {\partial h} {\partial t} + \operatorname{div} (h \bm u) = 0
\\[0.5 cm]
& \frac {\partial p} {\partial t} + \operatorname{div} (p \bm u) + p \operatorname{div} \bm u = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (h \bm u) + \operatorname{div} (h \bm u \otimes \bm u) + \operatorname{grad} p = 0
\end{aligned}
$$

'''