import streamlit as st


r"""
# Модельная задача

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
	&\frac{\partial \varrho}{\partial t} + \div(\varrho \bm u) = 0 ,
	\quad &x& \in \Omega ,
	\quad &0& < t \leq T
	\\[0.5 cm]
	&\frac{\partial }{\partial t} (\varrho \bm u) + \div(\varrho \bm u \otimes \bm u) + \grad p = 0 , 
	\quad &x& \in \Omega ,
	\quad &0& < t \leq T
	\\[0.5 cm]
	&\varrho(x, 0) = \varrho^0(x) ,
	\quad &x& \in \Omega
	\\[0.5 cm]
	&\bm u (x, 0) = \bm u^0(x) ,
	\quad &x& \in \Omega
	\\[1 cm]
	&\Omega = \{x \ | \ x = (x_1, x_2), \quad -5 < x_1 < 5, \quad -5 < x_2 < 5\}
	\\[0.5 cm]
	& a = 1
	\\[0.5 cm]
	& \gamma = 1.4
	\\[0.5 cm]
	&\bm u^0(x) = 0
	\\[0.5 cm]
	&\varrho^0(x) = 1 + \alpha \operatorname{exp}(-\beta |x|^2)
	\\[0.5 cm]
	&\alpha = 2
	\\[0.5 cm]
	&\beta = 20
\end{aligned}
$$
"""
