import streamlit as st


r"""
# Дифференциально-операторная формулировка

Для удобства рассмотрения введем операторы адвективного (конвективного) переноса в 
системе уравнений Эйлера. Оператор адвекции $\mathcal{A} = \mathcal{A}(\bm u)$
записывается в дивергентной форме следующим образом:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \mathcal{A}(\bm u) \varphi  =  \div  (\bm u \, \varphi ) . 
\end{aligned}
$$

При условии, что для скорости $\bm u$ выполнено граничное условие (3), имеет место

$$
\begin{aligned}
 (\mathcal{A}(\bm u) \varphi, 1) = 0 .
\end{aligned}
$$


Уравнение неразрывности (1) запишется в виде дифференциально-операторного уравнения

$$
\begin{aligned}
 \frac{d \varrho }{d t} + \mathcal{A}(\bm u) \varrho  = 0,
 \quad 0 < t \leq T, 
\end{aligned}
$$

при использовании обозначений $\varrho(t) = \varrho (x,t)$. 
Аналогично, уравнение (2) записывается в виде

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{d  }{d t} (\varrho \bm u) + \mathcal{A}(\bm u) (\varrho \bm u) + \grad p(\varrho) = 0,
 \quad 0 < t \leq T.  
\end{aligned} 
$$

Рассматривается задача Коши для системы уравнений, состоящией из двух предыдущих уравнений, при заданной зависимости
$p(\varrho)$, когда начальные условия имеют вид (4)

$$
\begin{aligned}
 \varrho(0) = \varrho^0,
 \quad \bm u (0) = \bm u^0 .
\end{aligned} 
$$

Для рассматриваемой задачи ключевое значение имеет свойство оператора адвекции, записанного в дивергентном виде.
"""
