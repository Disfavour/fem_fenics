import streamlit as st


r"""
# Вариационная постановка

Простой подход к решению зависящих от времени краевых задач методом конечных
элементов заключается в том, чтобы сначала дискретизировать производную по времени с помощью конечно-разностной
аппроксимации, что приводит к последовательности стационарных задач, а затем преобразовать
каждую стационарную задачу в вариационную формулировку.

Пусть верхний индекс $n$ обозначает величину в момент времени $t_n$, где $n$ - целое число, подсчитывающее уровни
времени. Например, $u^n$ означает $u$ на временном уровне $n$. Дискретизация с конечной разностью
во времени сначала состоит из выборки краевой задачи на некотором временном уровне $t_{n+1}$:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

{\left( \frac {\partial u} {\partial t} \right)}^{n+1} = \div \grad u^{n+1} + f^{n+1}.
$$

Производная по времени может быть аппроксимирована:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

{\left( \frac {\partial u} {\partial t} \right)}^{n+1} = \frac {u^{n+1} - u^{n}} {\Delta t},
$$

где $\Delta t$ - параметр дискретизации по времени. Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\frac {u^{n+1} - u^{n}} {\Delta t} = \div \grad u^{n+1} + f^{n+1}.
$$

Это дискретная по времени версия уравнения теплопроводности, так называемая обратная (неявная) дискретизация Эйлера.

Преобразуем таким образом, чтобы левая часть содержала компоненты с
неизвестным $u^{n+1}$, а правая часть содержала только вычисленные компоненты:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

u^{n+1} - \Delta t \div \grad u^{n+1} = u^n + \Delta t f^{n+1}
$$

Зная $u_0$ можно вычислить $u_1$ и так далее.

Необходимо получить вариационные формулировки для уравнений:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

u^{n+1} - \Delta t \div \grad u^{n+1} = u^n + \Delta t f^{n+1}.
$$

Умножаем обе части дифференциального уравнения на тестовую функцию $v$ и интегрируем по области:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\int \limits_\Omega (u^{n+1}v - \Delta t \div \grad u^{n+1}v) \ dx
= \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
$$

Рассмотрим следующую формулу векторного анализа:

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

После применения формулы к компоненту левой части уравнения получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega u^{n+1}v \ dx
- \int \limits_\Omega \Delta t \div \grad (u^{n+1}v) \ dx
+ \int \limits_\Omega \Delta t \grad u^{n+1} \grad v \ dx
= \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
\end{aligned}
$$

Рассмотрим теорему о дивергенции:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega \operatorname{div} \bold{f} \ dx = \int \limits_{\partial \Omega} \bold{f} \cdot n \ ds.
\end{aligned}
$$

Применим её к соответствующему компоненту из левой части уравнения:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega u^{n+1}v \ dx
- \Delta t \int \limits_{\partial \Omega} \grad (u^{n+1}v) \cdot n \ ds
+ \Delta t \int \limits_\Omega \grad u^{n+1} \grad v \ dx
= \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
\end{aligned}
$$

На границе, где задано условие Дирихле, можно выбирать тестовую функцию так, чтобы обнулить соответствующий интеграл:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\int \limits_\Omega u^{n+1}v \ dx
+ \Delta t \int \limits_\Omega \grad u^{n+1} \grad v \ dx
= \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
\end{aligned}
$$

Можно получить выражение в виде билинейной и линейной формах:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&a(u^{n+1}, \ v) = L_{n+1}(v),
\\[0.5 cm]
&a(u^{n+1}, \ v) = \int \limits_\Omega u^{n+1}v \ dx + \Delta t \int \limits_\Omega \grad u^{n+1} \grad v \ dx,
\\[0.5 cm]
&L_{n+1}(v) = \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
\end{aligned}
$$

Также существует альтернативная форма записи:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&F_{n+1}(u^{n+1}, \ v) = 0,
\\[0.5 cm]
&F_{n+1}(u^{n+1}, \ v) = \int \limits_\Omega u^{n+1}v \ dx + \Delta t \int \limits_\Omega \grad u^{n+1} \grad v \ dx
- \int \limits_\Omega (u^nv + \Delta t f^{n+1}v) \ dx.
\end{aligned}
$$

В дополнение к вариационной задаче, которая должна решаться на каждом временном шаге, нам
также необходимо аппроксимировать начальное условие. Это уравнение также может
быть преобразовано в вариационную задачу аналогичным образом:

$$
\begin{aligned}
&u^0 = u_0,
\\[0.5 cm]
&\int \limits_\Omega u^0 v \ dx = \int \limits_\Omega u_0 v \ dx,
\\[1.0 cm]
&a_0(u^0, \ v) = L_0(v),
\\[0.5 cm]
&a_0(u^0, \ v) = \int \limits_\Omega u^0 v \ dx,
\\[0.5 cm]
&L_0(v) = \int \limits_\Omega u_0 v \ dx.
\end{aligned}
$$

Таким образом, нам необходимо на каждом шаге по времени находить:

- $u^0$ из $a_0(u^0, \ v) = L_0(v)$;
- $u^{n+1}$ из $a(u^{n+1}, \ v) = L_{n+1}(v)$.


"""
