import streamlit as st


r"""
# Математическая модель

Уравнение неразрывности в ограниченной области $\Omega$ имеет вид

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
	\frac{\partial \varrho}{\partial t} + \div(\varrho \bm u) = 0 ,
	\quad x \in \Omega ,
	\quad 0 < t \leq T
\end{equation}
$$

где $\varrho(x, t) > 0$ - плотность, а $\bm u (x, t)$ - скорость.

Уравнение движения в консервативном виде, когда

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{equation}
	\frac{\partial }{\partial t} (\varrho \bm u) + \div(\varrho \bm u \otimes \bm u) + \grad p = 0 , 
	\quad x \in \Omega ,
	\quad 0 < t \leq T
\end{equation}
$$

где $p(x, t)$ - давление.

Жидкость предполагается баротропной, так что предполагается
известной зависимость давления от плотности: $p = p(\varrho)$, ${\displaystyle \frac{d p}{d \varrho} > 0}$.

Границы считаются твердыми. В силу этого
имеем граничное условие непротекания

$$
\begin{equation}
	(\bm u \cdot \bm n) = 0, 
	\quad x \in \partial \Omega 
\end{equation}
$$

Задаются также начальные условия для плотности и скорости:

$$
\begin{equation}
	\varrho(x, 0) = \varrho^0(x) ,
	\quad \bm u (x, 0) = \bm u^0(x) ,
	\quad x \in \Omega
\end{equation}
$$

Начально-краевая задача (1)-(4) описывает нестационарные течения 
идеальной баротропной жидкости.
"""
