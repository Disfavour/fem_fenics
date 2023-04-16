import streamlit as st


r"""
# Векторный анализ

- $\varphi,\ \psi, f_i$ - скалярные функции;
- $y$ - точка;
- $\bold {g}, \bold {f}$ - векторные поля;
- $a,\ b$ - вещественные числа.


Векторное поле:

$$
\bold{f} = \left( f_1, f_2, \dots, f_n \right).
$$



## Градиент

$$
\operatorname{grad} \varphi = \lim_{V(\Omega) \to 0} \frac {1} {V(\Omega)} \left( \int \limits_{\partial \Omega} \varphi \ ds
\right)
$$

Определение в координатах:

$$
\operatorname{grad} \varphi = \left(
\frac {\partial \varphi} {\partial x_1}, \dots, \frac {\partial \varphi} {\partial x_n} \right).
$$


## Дивергенция

$$
\operatorname{div} \bold{f} = \lim_{V(\Omega) \to 0} \frac {1} {V(\Omega)} \left( \int \limits_{\partial \Omega} \bold{f}
\cdot n \ ds
\right)
$$

Определение в координатах:

$$
\operatorname{div} \bold{f} = \frac {\partial f_1} {\partial x_1}
+ \frac {\partial f_2} {\partial x_2} + \dots + \frac {\partial f_n} {\partial x_n}.
$$


## Ротор

$$
\operatorname{rot} \bold{f} = \lim_{V(\Omega) \to 0} \frac {1} {V(\Omega)} \left( \int \limits_{\partial \Omega} n \times
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
\nabla = \left( \frac {\partial} {\partial x_1}, \dots, \frac {\partial} {\partial x_n} \right)
\\[0.5 cm]
\operatorname{grad} \varphi = \nabla \varphi
\\[0.5 cm]
\operatorname{div} \bold{f} = \nabla \cdot \bold{f}
\\[0.5 cm]
\operatorname{rot} \bold{f} = \nabla \times \bold{f}
$$


## Лапласиан

$$
\Delta f = \operatorname{div} \operatorname{grad} f = \nabla^2 f = \nabla \cdot \nabla f
= \sum \limits_{i=1}^{n} \frac {\partial^2 f} {\partial x_i^2}
$$


## Теорема о дивергенции (Остроградского-Гаусса)

$$
\int \limits_\Omega \operatorname{div} \bold{f} \ d\bm{x} = \int \limits_{\partial \Omega} \bold{f} \cdot n \ ds,
$$

где $n$ - внешняя нормаль.


### Нормальная производная

$$
\frac {\partial u} {\partial n} = \operatorname{grad} u \cdot n
$$
"""
