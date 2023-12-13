import streamlit as st

r'''
## Координатная форма

Получим координатную форму уравнения мелкой воды из дивергентного
$$
\begin{aligned}
& \frac {\partial \varrho} {\partial t} + \operatorname{div} (\varrho \bm u) = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (\varrho \bm u) + \operatorname{div} (\varrho \bm u \otimes \bm u)
+ \operatorname{grad} p = 0
\end{aligned}
$$

Уравнение неразрывности в координатном виде
$$
\frac {\partial \varrho} {\partial t}
+ \frac {\partial} {\partial x_1} (\varrho u)
+ \frac {\partial} {\partial x_2} (\varrho v) = 0
$$

В уравнении движения раскроем векторы и внешнее произведение
$$
\frac {\partial} {\partial t}
\begin{pmatrix}
   \varrho u \\[0.3cm]
   \varrho v
\end{pmatrix}
+ \operatorname{div}
\begin{pmatrix}
   \varrho u^2 & \varrho uv \\[0.3cm]
   \varrho uv & \varrho v^2
\end{pmatrix}
+
\begin{pmatrix}
   \displaystyle \frac {\partial p} {\partial x_1} \\[0.3cm]
   \displaystyle \frac {\partial p} {\partial x_2}
\end{pmatrix}
= 0
$$

Раскроем дивергенцию
$$
\begin{pmatrix}
   \displaystyle \frac {\partial} {\partial t} (\varrho u) \\[0.3cm]
   \displaystyle \frac {\partial} {\partial t} (\varrho v)
\end{pmatrix}
+
\begin{pmatrix}
   \displaystyle \frac {\partial} {\partial x_1} (\varrho u^2) + \displaystyle \frac {\partial} {\partial x_2} (\varrho uv) \\[0.3cm]
   \displaystyle \frac {\partial} {\partial x_1} (\varrho uv) + \displaystyle \frac {\partial} {\partial x_2} (\varrho v^2)
\end{pmatrix}
+
\begin{pmatrix}
   \displaystyle \frac {\partial p} {\partial x_1} \\[0.3cm]
   \displaystyle \frac {\partial p} {\partial x_2}
\end{pmatrix}
= 0
$$

Можем записать уравнения отдельно
$$
\begin{aligned}
& \frac {\partial} {\partial t} (\varrho u)
+ \frac {\partial} {\partial x_1} (\varrho u^2)
+ \frac {\partial} {\partial x_2} (\varrho uv)
+ \frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \frac {\partial} {\partial t} (\varrho v)
+ \frac {\partial} {\partial x_1} (\varrho uv)
+ \frac {\partial} {\partial x_2} (\varrho v^2)
+ \frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$

Дифференцируем
$$
\begin{aligned}
& \left( \varrho \frac {\partial u} {\partial t} + u \frac {\partial \varrho} {\partial t}\right)
+ \left( u \frac {\partial} {\partial x_1} (\varrho u) + \varrho u \frac {\partial u} {\partial x_1}\right)
+ \left( u \frac {\partial} {\partial x_2} (\varrho v) + \varrho v \frac {\partial u} {\partial x_2}\right)
+ \frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \left( \varrho \frac {\partial v} {\partial t} + v \frac {\partial \varrho} {\partial t}\right)
+ \left( v \frac {\partial} {\partial x_1} (\varrho u) + \varrho u \frac {\partial v} {\partial x_1}\right)
+ \left( v \frac {\partial} {\partial x_2} (\varrho v) + \varrho v \frac {\partial v} {\partial x_2}\right)
+ \frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$

Используем уравнение неразрывности
$$
\begin{aligned}
& \varrho \frac {\partial u} {\partial t} + u \left( -\frac {\partial} {\partial x_1} (\varrho u) - \frac {\partial} {\partial x_2} (\varrho v) \right)
+ u \frac {\partial} {\partial x_1} (\varrho u) + \varrho u \frac {\partial u} {\partial x_1}
+ u \frac {\partial} {\partial x_2} (\varrho v) + \varrho v \frac {\partial u} {\partial x_2}
+ \frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \varrho \frac {\partial v} {\partial t} + v \left( -\frac {\partial} {\partial x_1} (\varrho u) - \frac {\partial} {\partial x_2} (\varrho v) \right)
+ v \frac {\partial} {\partial x_1} (\varrho u) + \varrho u \frac {\partial v} {\partial x_1}
+ v \frac {\partial} {\partial x_2} (\varrho v) + \varrho v \frac {\partial v} {\partial x_2}
+ \frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$

Сокращаем слагаемые и делим на $\varrho$
$$
\begin{aligned}
& \frac {\partial u} {\partial t}
+ u \frac {\partial u} {\partial x_1}
+ v \frac {\partial u} {\partial x_2}
+ \frac {1} {\varrho} \frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \frac {\partial v} {\partial t}
+ u \frac {\partial v} {\partial x_1}
+ v \frac {\partial v} {\partial x_2}
+ \frac {1} {\varrho} \frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$

Таким образом, уравнение мелкой воды в координатной форме принимает вид
$$
\begin{aligned}
& \frac {\partial \varrho} {\partial t}
+ \frac {\partial} {\partial x_1} (\varrho u)
+ \frac {\partial} {\partial x_2} (\varrho v) = 0
\\[0.5 cm]
& \frac {\partial u} {\partial t}
+ u \frac {\partial u} {\partial x_1}
+ v \frac {\partial u} {\partial x_2}
+ \frac {1} {\varrho} \frac {\partial p} {\partial x_1} = 0
\\[0.5 cm]
& \frac {\partial v} {\partial t}
+ u \frac {\partial v} {\partial x_1}
+ v \frac {\partial v} {\partial x_2}
+ \frac {1} {\varrho} \frac {\partial p} {\partial x_2} = 0
\end{aligned}
$$
'''
