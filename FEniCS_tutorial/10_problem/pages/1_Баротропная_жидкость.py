import streamlit as st


r"""
# Баротропная жидкость

Уравнение неразрывности в ограниченной области $\Omega$ имеет вид

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
 \frac{\partial \varrho}{\partial t} + \div(\varrho \bm u) = 0 ,
 \quad x \in \Omega ,
 \quad 0 < t \leq T ,
\end{equation}
$$

где $\varrho(x, t) > 0$ - плотность, а $\bm u (x, t)$ - скорость.
Уравнение движения запишем в консервативном виде, когда

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
 \frac{\partial }{\partial t} (\varrho \bm u) + \div(\varrho \bm u \otimes \bm u) + \grad p = 0 , 
 \quad x \in \Omega ,
 \quad 0 < t \leq T ,
\end{equation}
$$

где $p(x, t)$ - давление. Жидкость предполагается баротропной, так что предполагается
известной зависимость давления от плотности: $p = p(\varrho)$, ${\displaystyle \frac{d p}{d \varrho} > 0}$.

Границы считаются твердыми. В силу этого
имеем граничное условие непротекания

$$
\begin{equation}
 (\bm u \cdot \bm n) = 0, 
 \quad x \in \partial \Omega . 
\end{equation}
$$

Задаются также начальные условия для плотности и скорости:

$$
\begin{equation}
 \varrho(x, 0) = \varrho^0(x) ,
 \quad \bm u (x, 0) = \bm u^0(x) ,
 \quad x \in \Omega .
\end{equation}
$$

Начально - краевая задача (1)-(4) описывает нестационарные течения 
идеальной баротропной жидкости.

Непосредственным интегрированием уравнения неразрывности (1) по области $\Omega$ с учетом граничного условия (3) можно 
получить закон сохранения массы:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx + \int \limits_\Omega \div(\varrho \bm u) \ dx = 0
\end{aligned}
$$

После применения теоремы о дивергенции получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx
+ \int \limits_{\partial \Omega} \varrho \bm u \cdot n \ dx = 0
\end{aligned}
$$

С учетом граничных условий получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial \varrho} {\partial t} \ dx = 0
\\[0.5 cm]
&\frac {\partial} {\partial t} \int \limits_\Omega \varrho \ dx = 0
\end{aligned}
$$

Таким образом, масса не зависит от времени:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \varrho(x, t) \ dx = \int \limits_\Omega \varrho(x, 0) \ dx = \int \limits_\Omega \varrho^0(x) \ dx
\end{aligned}
$$

$$
\begin{equation}
 m (t) = m(0) ,
 \quad m(t) = \int_{\Omega}  \varrho(x, t) d x .
\end{equation}
$$

В гильбертовом пространстве $L_2(\Omega)$ мы определяем скалярное произведение и норму стандартным способом:

$$
  (w,u) = \int_{\Omega} w({x}) u({x}) d{x},
  \quad \|w\| = (w,w)^{1/2} .
$$

Аналогично определяется пространство векторных функций $\bm L_2(\Omega)$.
При неотрицательности плотности закон сохранения массы можно записать в виде

$$
 \|\varrho^{1/2}\| = \|\varrho_0^{1/2}\| .
$$

Это соотношение можно рассматривать как априорную оценку в $L_2(\Omega)$ для $\varrho^{1/2}$. 

Уравнение (2) напрямую выражает закон созранения импульса. 
Интегрируя это уравнение по $\Omega$, получим

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho u) \ dx
+ \int \limits_\Omega \div (\varrho u \otimes u) \ dx
+ \int \limits_\Omega \grad p \ dx = 0
\end{aligned}
$$

После применения теорем о дивергенции и градиенте получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho u) \ dx
+ \int \limits_{\partial \Omega} \varrho u \otimes u \cdot n \ dx
+ \int \limits_{\partial \Omega} p \bm n \ dx = 0
\end{aligned}
$$

С учетом граничных условий получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\int \limits_\Omega \frac {\partial} {\partial t} (\varrho u) \ dx
+ \int \limits_{\partial \Omega} p \bm n \ dx = 0\end{aligned}
$$

Таким образом, мы имеем:

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
\begin{equation}
 \bm I (t) = \bm I(0) - \int_{\partial \Omega} p(\varrho) \bm n d x ,
 \quad \bm I(t) =  \int_{\Omega} \varrho \bm u d x .
\end{equation}
$$

Домножая на $\bm u$ и принимая во внимание уравнение (1), запишем уравнение (2) в виде

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\bm u \frac {\partial} {\partial t} (\varrho \bm u) + u \div (\varrho \bm u \otimes \bm u) + \bm u \grad p = 0
\\[0.5 cm]
&\bm u \left( \varrho \frac {\partial \bm u} {\partial t} + \bm u \frac {\partial \varrho} {\partial t} \right)
+ \bm u \div (\varrho \bm u \otimes \bm u) + (\div(p \bm u) - p \div \bm u) = 0
\end{aligned}
$$

Учтем:

$$
\frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2)
= \varrho \bm u \frac {\partial \bm u} {\partial t} + \frac {1} {2} \bm u^2 \frac {\partial \varrho} {\partial t}
$$

Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2)
+ \frac {1} {2} \bm u^2 \frac {\partial \varrho} {\partial t}
+ \bm u \div (\varrho \bm u \otimes \bm u) + \div(p \bm u) - p \div \bm u = 0
\end{aligned}
$$

С учетом уравнения непрерывности:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2)
- \frac {1} {2} \bm u^2 \div (\varrho \bm u)
+ \bm u \div (\varrho \bm u \otimes \bm u) + \div(p \bm u) - p \div \bm u = 0
\end{aligned}
$$

Учтем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {1} {2} \div (\varrho \bm u^3)
= \frac {1} {2} ( \varrho \bm u \cdot \grad \bm u^2 + \bm u^2 \div (\varrho \bm u))
\end{aligned}
$$

Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2)
- \frac {1} {2} \div (\varrho \bm u^3) + \frac {1} {2} \varrho \bm u \cdot \grad \bm u^2
+ \bm u \div (\varrho \bm u \otimes \bm u) + \div(p \bm u) - p \div \bm u = 0
\end{aligned}
$$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
&\frac {1} {2} \frac {\partial} {\partial t} (\varrho \bm u^2)
+ \frac {1} {2} \div (\varrho \bm u^3)
+ \div(p \bm u) - p \div \bm u = 0
\end{aligned}
$$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

 \frac{1}{2} \frac{\partial }{\partial t} (\varrho |\bm u |^2) +
 \frac{1}{2} \div (\varrho |\bm u |^2 \bm u) + \div (p(\varrho) \bm u) - p \div \bm u = 0 
$$

Интегрируем по области:

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

После применения теоремы о дивергенции получаем:

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

С учетом граничных условий:

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

Таким образом, интегрирование по области $\Omega$ с учетом (3) дает

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
 \frac{1}{2} \frac{d }{d t} \int_{\Omega} \varrho |\bm u |^2 d x - 
 \int_{\Omega} p(\varrho) \div \bm u \, d x = 0 .
\end{equation}
$$

Определим потенциал давления $\varPi(\varrho)$ из уравнения

$$
\begin{equation}
 \varrho \frac{d \varPi}{d \varrho} - \varPi(\varrho) = p(\varrho) .
\end{equation} 
$$

В частности, для идеальной жидкости имеем

$$
\begin{equation}
 p(\varrho)  = a\varrho^\gamma, 
 \quad \varPi (\varrho) = a \frac{\varrho^\gamma}{\gamma -1} ,
 \quad a = \operatorname{const} > 0,
 \quad \gamma > 1 . 
\end{equation} 
$$

Из уравнения неразрывности (1) имеем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
 \frac{\partial \varPi }{\partial t} + \div(\varPi  \bm u) + p(\varrho) \div \bm u = 0 ,
 \quad 0 < t \leq T . 
\end{equation} 
$$

Интегрирование (10) дает

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

 \frac{d }{d t} \int_{\Omega} \varPi \, d x +
 \int_{\Omega} p(\varrho) \div \bm u \, d x = 0 .
$$

Складывая это равенство с (7) получим

$$
\begin{equation}
 \frac{d }{d t}  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 
 + \varPi(\varrho) \right ) d x = 0 .
\end{equation}
$$

Приходим к закону сохранения полной механической энергии

$$
\begin{equation}
 E(t) = E(0),
 \quad E(t) =  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 
 + \varPi(\varrho) \right ) d x .
\end{equation}
$$

Равенства (5), (6) и (12) являются основными законами сохранения для 
задачи (1)-(4). 
"""