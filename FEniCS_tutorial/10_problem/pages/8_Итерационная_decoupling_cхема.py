import streamlit as st


r"""
# Итерационная decoupling cхема

На основе линеаризованной схемы 
можно строить итерационный алгоритм для вычислительной реализации чисто неявной схемы.
Приближенное решение для $\varrho_{n+1}$, $\bm u_{n+1}$ на $k$-oй итерации обозначим
$\varrho_{n+1}^{k}$, $\bm u_{n+1}^{k}$, причем используется начальное приближение 
по предыдущему слою по времени:

$$
\begin{aligned}
 \varrho_{n+1}^0 = \varrho_{n},
 \quad \bm u_{n+1}^0 = \bm u_{n} .
\end{aligned} 
$$

Будем считать, что новое приближение  на новом слое по времени находится, когда делается $K$ итераций. 
Аналогично используется система уравнений

$$
\begin{aligned}
 \frac{\varrho_{n+1}^{k+1} - \varrho_{n} }{\tau} + A(\bm u_{n+1}^k) \varrho_{n+1}^{k+1}  = 0,
\end{aligned} 
$$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{\varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1} - \varrho_{n} \bm u_{n}  }{\tau } + 
 A(\bm u_{n+1}^{k}) (\varrho_{n+1}^{k+1} \bm u_{n+1}^{k+1}) + \grad p(\varrho_{n+1}^{k+1}) = 0, 
 \quad k = 0, 1, ..., K-1 ,
\end{aligned} 
$$

причем

$$
\begin{aligned}
 \varrho_{n+1} = \varrho_{n+1}^K,
 \quad \bm u_{n+1} = \bm u_{n+1}^K ,
 \quad n = 0, 1, ..., N-1 .
\end{aligned} 
$$

Тем самым, на каждой итерации сначала решается линейная задача для плотности, а потом - линейная
задача для скорости.
"""
