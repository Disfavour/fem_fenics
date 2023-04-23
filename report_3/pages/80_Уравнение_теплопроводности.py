import streamlit as st


r"""
# Уравнение теплопроводности

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
-&\div \grad u +\frac {\partial u} {\partial t} = f,         \quad &x& \in \Omega, \ &t \in (0, T]
\\[0.5 cm]
&u = u_D,                                                    \quad &x& \in \partial \Omega, \ &t \in (0, T]
\\[0.5 cm]
&u = u_0, \quad &t& = 0
\end{aligned}
$$

Здесь $u$ изменяется в зависимости от пространства и времени ( $u = u(x_1, x_2, \dots, t)$). Исходная функция $f$ и
граничные значения $u_D$ также могут изменяться в зависимости от пространства и времени. Начальное условие $u_0$
является функцией только пространства.

Простой подход к решению зависящих от времени краевых задач методом конечных
элементов заключается в том, чтобы сначала дискретизировать производную по времени с помощью конечно-разностной
аппроксимации, что приводит к последовательности стационарных задач, а затем преобразовать
каждую стационарную задачу в вариационную формулировку.

Для приближенного решения будет использовать $y$, для которого выполнено:

$$
y \approx u
$$
"""
