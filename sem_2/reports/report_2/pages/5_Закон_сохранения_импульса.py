import streamlit as st

r'''
## Закон сохранения импульса

Уравнение движения выражает закон сохранения импульса

Интегрируем уравнение по $\Omega$
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho \bm u) \ dx
+ \int \limits_\Omega \div (\varrho \bm u \otimes \bm u) \ dx
+ \int \limits_\Omega \grad p \ dx = 0
\end{aligned}
$$

Применяем теоремы о дивергенции и градиенте
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho \bm u) \ dx
+ \int \limits_{\partial \Omega} \varrho \bm u \otimes \bm u \cdot n \ dx
+ \int \limits_{\partial \Omega} p n \ dx = 0
\end{aligned}
$$

Учтем граничное условие
$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho \bm u) \ dx
+ \int \limits_{\partial \Omega} p n \ dx = 0\end{aligned}
$$

Таким образом, получаем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&I = \int \limits_\Omega (\varrho \bm u) \ dx
\\[0.5 cm]
& \int \limits_0^T \frac {\partial I} {\partial t} dt = I(t) - I(0)
\end{aligned}
$$

Тем самым имеет место 

$$
\begin{aligned}
&\bm I (t) = \bm I(0) - \int_{\partial \Omega} p(\varrho) n d x
\\[0.5 cm]
&\bm I(t) =  \int_{\Omega} \varrho \bm u d x
\end{aligned}

$$
'''
