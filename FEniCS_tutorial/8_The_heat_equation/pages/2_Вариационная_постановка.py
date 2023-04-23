import streamlit as st


r"""
# Вариационная постановка

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

-\div \grad y^{n+1} + \frac {y^{n+1} - y^{n}} {\tau} \approx f^{n+1}
$$

Преобразуем таким образом, чтобы левая часть содержала компоненты с
неизвестным $u^{n+1}$, а правая часть содержала только вычисленные компоненты:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

- \tau \div \grad y^{n+1} + y^{n+1}  = y^n + \tau f^{n+1}
$$

Таким образом, $y^1$ может быть вычислен по $y^0$ и так далее.

Умножаем обе части дифференциального уравнения на тестовую функцию $v$ и интегрируем по области:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

- \int \limits_\Omega \tau \div \grad y^{n+1} v \ dx
+ \int \limits_\Omega y^{n+1} v \ dx
= \int \limits_\Omega y^n v \ dx
+ \int \limits_\Omega \tau f^{n+1} v \ dx
$$

Применим следующую формулу векторного анализа:

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

Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- \int \limits_\Omega \tau \div \grad (y^{n+1} v) \ dx
+ \int \limits_\Omega \tau \grad y^{n+1} \grad v \ dx
+ \int \limits_\Omega y^{n+1} v \ dx
= \int \limits_\Omega y^n v \ dx
+ \int \limits_\Omega \tau f^{n+1} v \ dx
\end{aligned}
$$

Применим теорему о дивергенции:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega \operatorname{div} \bold{f} \ dx = \int \limits_{\partial \Omega} \bold{f} \cdot n \ ds
\end{aligned}
$$

Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
- \int \limits_{\partial \Omega} \tau \grad (y^{n+1} v) \cdot n \ dx
+ \int \limits_\Omega \tau \grad y^{n+1} \grad v \ dx
+ \int \limits_\Omega y^{n+1} v \ dx
= \int \limits_\Omega y^n v \ dx
+ \int \limits_\Omega \tau f^{n+1} v \ dx
\end{aligned}
$$

На границе, где задано условие Дирихле, можно выбирать тестовую функцию так, чтобы обнулить соответствующий интеграл:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega \tau \grad y^{n+1} \grad v \ dx
+ \int \limits_\Omega y^{n+1} v \ dx
= \int \limits_\Omega y^n v \ dx
+ \int \limits_\Omega \tau f^{n+1} v \ dx
\end{aligned}
$$

Выражение в виде билинейной и линейной формах:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&a(y^{n+1}, \ v) = \int \limits_\Omega \tau \grad y^{n+1} \grad v \ dx
+ \int \limits_\Omega y^{n+1} v \ dx
\\[0.5 cm]
&L(v) = \int \limits_\Omega y^n v \ dx
+ \int \limits_\Omega \tau f^{n+1} v \ dx
\end{aligned}
$$

В дополнение к вариационной задаче, которая должна решаться на каждом временном шаге, нам
также необходимо аппроксимировать начальное условие. Это уравнение также может
быть преобразовано в вариационную задачу аналогичным образом:

$$
\begin{aligned}
&y^0 \approx u^0 = u_0
\\[0.5 cm]
&\int \limits_\Omega y^0 v \ dx = \int \limits_\Omega u_0 v \ dx
\end{aligned}
$$
"""
