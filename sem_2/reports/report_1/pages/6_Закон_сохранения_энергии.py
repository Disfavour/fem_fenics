import streamlit as st

r'''
## Закон сохранения энергии

Домножая на $\bm u$ и принимая во внимание уравнение неразрывности, запишем уравнение движения в виде
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

 \frac{1}{2} \frac{\partial }{\partial t} (\varrho |\bm u |^2) +
 \frac{1}{2} \div (\varrho |\bm u |^2 \bm u) + \div (p(\varrho) \bm u) - p \div \bm u = 0 
$$

Интегрируем по области
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
& \int \limits_\Omega \frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2) \ dx
+ \int \limits_\Omega \frac {1} {2} \div (\varrho \bm u^3) \ dx
+ \int \limits_\Omega \div(p \bm u) \ dx - \int \limits_\Omega p \div \bm u \ dx = 0
\end{aligned}
$$

Применяем теорему о дивергенции
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
& \int \limits_\Omega \frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2) \ dx
+ \int \limits_{\partial \Omega} \frac {1} {2} \varrho \bm u^3 \cdot n \ dx
+ \int \limits_{\partial \Omega} p \bm u \cdot n \ dx - \int \limits_\Omega p \div \bm u \ dx = 0
\end{aligned}
$$

Учтем граничное условие

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
& \int \limits_\Omega \frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2) \ dx
- \int \limits_\Omega p \div \bm u \ dx = 0
\\[0.5 cm]
& \frac {1} {2} \frac {\partial} {\partial t} \int \limits_\Omega \varrho \bm u^2 \ dx
- \int \limits_\Omega p \div \bm u \ dx = 0
\end{aligned}
$$

Таким образом, получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{1}{2} \frac{d }{d t} \int_{\Omega} \varrho |\bm u |^2 d x - 
 \int_{\Omega} p(\varrho) \div \bm u \, d x = 0
\end{aligned}
$$

Определим потенциал давления $\varPi(\varrho)$ из уравнения

$$
\begin{aligned}
\varrho \frac{d \varPi}{d \varrho} - \varPi(\varrho) = p(\varrho) 
\end{aligned}
$$

В частности, для идеальной жидкости имеем

$$
\begin{aligned}
p(\varrho)  = a\varrho^\gamma, 
\quad \varPi (\varrho) = a \frac{\varrho^\gamma}{\gamma -1} ,
\quad a = \operatorname{const} > 0,
\quad \gamma > 1
\end{aligned}
$$

Из уравнения неразрывности имеем
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{\partial \varPi }{\partial t} + \div(\varPi  \bm u) + p(\varrho) \div \bm u = 0 ,
 \quad 0 < t \leq T
\end{aligned}
$$

Интегрируем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

 \frac{d }{d t} \int_{\Omega} \varPi \, d x +
 \int_{\Omega} p(\varrho) \div \bm u \, d x = 0
$$

Складывая это равенство с предыдущим результатом получим

$$
\begin{aligned}
 \frac{d }{d t}  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 
 + \varPi(\varrho) \right ) d x = 0
\end{aligned}
$$

Приходим к закону сохранения полной механической энергии

$$
\begin{aligned}
 E(t) = E(0),
 \quad E(t) =  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 
 + \varPi(\varrho) \right ) d x
\end{aligned}
$$
'''
