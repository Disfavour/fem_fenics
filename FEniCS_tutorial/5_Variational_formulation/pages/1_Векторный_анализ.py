import streamlit as st


r"""
# Векторный анализ

Векторное поле:

$$
\bold{f} = \left( f_1, f_2, \dots, f_n \right).
$$

Пусть $\varphi,\ \psi$ - скалярная функция; $y$ - точка; $\bold {g}$ - векторное поле; $a,\ b$ - вещественные числа.


## Градиент

$$
\operatorname{grad} \varphi \big|_{y} = \lim_{V \to 0} \frac {1} {V} \left( \int \limits_{\partial \Omega} \varphi \ ds
\right)
$$

Определение в координатах:

$$
\operatorname{grad} \varphi = \left(
\frac {\partial \varphi} {\partial x_1}, \dots, \frac {\partial \varphi} {\partial x_n} \right).
$$


## Дивергенция

$$
\operatorname{div} \bold{f} \big|_{y} = \lim_{V \to 0} \frac {1} {V} \left( \int \limits_{\partial \Omega} \bold{f}
\cdot \bm n \ ds
\right)
$$

Определение в координатах:

$$
\operatorname{div} \bold{f} = \frac {\partial f_1} {\partial x_1}
+ \frac {\partial f_2} {\partial x_2} + \dots + \frac {\partial f_n} {\partial x_n}.
$$


## Ротор

$$
\operatorname{rot} \bold{f} \big|_{y} = \lim_{V \to 0} \frac {1} {V} \left( \int \limits_{\partial \Omega} \bm n \times
\bold{f} \ ds \right)
$$

В трёхмерной декартовой системе координат ротор вычисляется следующим образом:

$$
\operatorname{rot} \bold{f} =
\begin{vmatrix}
   \bold{i} & \bold{j} & \bold{k}
   \\[0.3 cm]
   \displaystyle \frac {\partial} {\partial x_1} & \displaystyle \frac {\partial} {\partial x_2}
   & \displaystyle \frac {\partial} {\partial x_3}
   \\[0.3 cm]
   f_1 & f_2 & f_3
\end{vmatrix}
= \left( \frac {\partial f_3} {\partial x_2} - \frac {\partial f_2} {\partial x_3};\ 
\frac {\partial f_1} {\partial x_3} - \frac {\partial f_3} {\partial x_1};\ 
\frac {\partial f_2} {\partial x_1} - \frac {\partial f_1} {\partial x_2} \right).
$$


## Набла

$$
\operatorname{grad} \varphi = \nabla \varphi
\\[0.5 cm]
\operatorname{div} \bold{f} = \nabla \cdot \bold{f}
\\[0.5 cm]
\operatorname{rot} \bold{f} = \nabla \times \bold{f}
\\[0.5 cm]
\nabla = \left( \frac {\partial} {\partial x_1}, \dots, \frac {\partial} {\partial x_n} \right)
$$


## Лапласиан

$$
\Delta f = \operatorname{div} \operatorname{grad} f = \nabla^2 f = \nabla \cdot \nabla f
= \sum \limits_{i=1}^{n} \frac {\partial^2 f} {\partial x_i^2}
$$


## Теорема о дивергенции (формула Гаусса-Остроградского)

$$
\int \limits_\Omega \operatorname{div} \bold{f} dx = \int \limits_{\partial \Omega} \bold{f} \cdot \bm{n} \ ds
$$


## Первая формула Грина

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}

\int \limits_\Omega (\psi \div \grad \varphi + \grad \psi \cdot \grad \varphi) \ dx
= \int \limits_{\partial \Omega} \psi (\grad \varphi \cdot \bm n) \ ds
= \int \limits_{\partial \Omega} \psi \frac {\partial \varphi} {\partial \bm{n}}\ ds
$$

где, $\displaystyle \frac {\partial} {\partial \bm{n}}$ - нормальная производная по направлению внешней нормали.
"""
