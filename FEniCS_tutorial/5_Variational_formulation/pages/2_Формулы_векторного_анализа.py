import streamlit as st


r"""
# Формулы векторного анализа

- $a$ - вещественное число.
- $\varphi,\ \psi$ - скалярные функции.
- $\bold{f},\ \bold{g}$ - векторные поля.

$$
\def \phi {\varphi}
\def \f {\bold{f}}
\def \g {\bold{g}}

\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
& \grad (a \psi + \phi) = a \grad \psi + \grad \phi
\\[0.5 cm]
& \div (a \f + \g) = a \div \f + \div \g
\\[0.5 cm]
& \rot (a \f + \g) = a \rot \f + \rot \g
\\[1.0 cm]
& \rot \grad \phi = 0
\\[0.5 cm]
& \div \rot \f = 0
\\[0.5 cm]
& \rot \rot \f = \grad \div \f - \div \grad \f
\\[1.0 cm]
& \div (\phi \f) = \f \cdot \grad \phi + \phi \div \f
\\[0.5 cm]
& \rot (\phi \f) = \grad \phi \times \f + \phi \rot \f
\\[0.5 cm]
& \grad (\f \cdot \g) = \g \grad \f + \f \grad \g + \f \times \rot \g + \g \times \rot \f
\\[0.5 cm]
& \div (\f \times \g) = \g \cdot \rot \f - \f \cdot \rot \g
\\[0.5 cm]
& \rot (\f \times \g) = \f \div \g - \g \div \f + \f \grad \g - \g \grad \f
\end{aligned}
$$
"""
