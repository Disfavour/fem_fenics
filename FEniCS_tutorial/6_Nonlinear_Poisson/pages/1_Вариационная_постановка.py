import streamlit as st


r"""
# Вариационная постановка

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (q(u) \grad u) = f 
\end{aligned}
$$

Умножаем обе части дифференциального уравнения на тестовую функцию $v$ и интегрируем по области

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_\Omega v \div (q(u) \grad u) \ d\bm{x} = \int \limits_\Omega f v \ d\bm{x}
\end{aligned}
$$

Рассмотрим следующую формулу векторного анализа

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\def \f {\bold{f}}

\begin{aligned}
& \div (\varphi \f) = \f \cdot \grad \varphi + \varphi \div \f
\\[0.5 cm]
& \varphi \div \f = \div (\varphi \f) - \f \cdot \grad \varphi
\end{aligned}
$$

После применения формулы к левой части уравнения получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_\Omega \div (v q(u) \grad u) \ d\bm{x}
+ \int \limits_\Omega q(u) \grad u \cdot \grad v \ d\bm{x}
= \int \limits_\Omega f v \ d\bm{x}
\end{aligned}
$$

Теорема о дивергенции

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega \operatorname{div} \bold{f} \ d\bm{x} = \int \limits_{\partial \Omega} \bold{f} \cdot n \ d\bm{s}
\end{aligned}
$$

После применения теоремы о дивергенции к левому слагаемому получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_{\partial \Omega} v q(u) \grad u \cdot n \ d\bm{x}
+ \int \limits_\Omega q(u) \grad u \cdot \grad v \ d\bm{x}
= \int \limits_\Omega f v \ d\bm{x}
\end{aligned}
$$

На границе, где задано условие Дирихле, выбираем тестовую функцию так, чтобы обнулить соответствующий интеграл

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega q(u) \grad u \cdot \grad v \ d\bm{x}
= \int \limits_\Omega f v \ d\bm{x}
\\[0.5 cm]
&\int \limits_\Omega q(u) \grad u \cdot \grad v - f v \ d\bm{x} = 0
\end{aligned}
$$

Представление для FEniCS нелинейной задачи

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&F(u, v) = 0
\\[0.5 cm]
&F(u, v) = \int \limits_\Omega q(u) \grad u \cdot \grad v - f v \ d\bm{x} = 0
\end{aligned}
$$
"""
