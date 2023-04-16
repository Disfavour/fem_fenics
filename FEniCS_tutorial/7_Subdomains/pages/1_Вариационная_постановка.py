import streamlit as st


r"""
# Вариационная постановка

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\div (\kappa(x) \grad u) = f 
\end{aligned}
$$

Умножаем обе части дифференциального уравнения на тестовую функцию $v$ и интегрируем по области

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_\Omega v \div (\kappa(x) \grad u) \ dx = \int \limits_\Omega f v \ dx.
\end{aligned}
$$

Рассмотрим следующую формулу векторного анализа

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}
\def \f {\bold{f}}

\begin{aligned}
& \div (\varphi \f) = \f \cdot \grad \varphi + \varphi \div \f,
\\[0.5 cm]
& \varphi \div \f = \div (\varphi \f) - \f \cdot \grad \varphi.
\end{aligned}
$$

После применения формулы к левой части уравнения получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_\Omega \div (v \kappa(x) \grad u) \ dx
+ \int \limits_\Omega \kappa(x) \grad u \cdot \grad v \ dx
= \int \limits_\Omega f v \ dx.
\end{aligned}
$$

Теорема о дивергенции

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega \operatorname{div} \bold{f} \ dx = \int \limits_{\partial \Omega} \bold{f} \cdot n \ ds.
\end{aligned}
$$

После применения теоремы о дивергенции к левому слагаемому получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- &\int \limits_{\partial \Omega} v \kappa(x) \grad u \cdot n \ dx
+ \int \limits_\Omega \kappa(x) \grad u \cdot \grad v \ dx
= \int \limits_\Omega f v \ dx.
\end{aligned}
$$

На границе, где задано условие Дирихле, выбираем тестовую функцию так, чтобы обнулить соответствующий интеграл

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int  \limits_\Omega \kappa(x) \grad u \cdot \grad v  \ dx
= \int \limits_\Omega f v \ dx.
\end{aligned}
$$

После преобразования получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&F(u, v) = \int \limits_\Omega ( \kappa(x) \grad u \cdot \grad v - f v ) \ dx = 0.
\end{aligned}
$$

Если $\kappa$ - константа, то:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&a(u, v) = \int \limits_\Omega \kappa(x) \grad u \cdot \grad v \ dx,
\\[0.5 cm]
&L(v) = \int \limits_\Omega f v  \ dx.
\end{aligned}
$$
"""
