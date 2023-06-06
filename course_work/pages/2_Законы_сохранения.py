import streamlit as st


r"""
# Законы сохранения

## Закон сохранения массы

Непосредственным интегрированием уравнения неразрывности по области $\Omega$ с учетом граничного условия можно
получить закон сохранения массы:

$$
\begin{aligned}
	m (t) = m(0) ,
	\quad m(t) = \int_{\Omega}  \varrho(x, t) d x
\end{aligned}
$$

## Закон сохранения импульса

Уравнение движения напрямую выражает закон сохранения импульса. 
Интегрируя это уравнение по $\Omega$, получим

$$
\int_{\Omega} \frac{\partial }{\partial t} (\varrho \bm u) d x + \int_{\partial \Omega} p(\varrho) \bm n d x = 0
$$

Тем самым имеет место

$$
\begin{aligned}
	\bm I (t) = \bm I(0) - \int_{\partial \Omega} p(\varrho) \bm n d x ,
	\quad \bm I(t) =  \int_{\Omega} \varrho \bm u d x
\end{aligned}
$$

## Закон сохранения полной механической энергии

$$
\begin{aligned}
	E(t) = E(0),
	\quad E(t) =  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 
	+ \Pi(\varrho) \right ) d x
\end{aligned}
$$

Потенциал давления $\Pi(\varrho)$ определяется из уравнения

$$
\begin{aligned}
	\varrho \frac{d \Pi}{d \varrho} - \Pi(\varrho) = p(\varrho)
\end{aligned} 
$$

где

$$
\begin{aligned}
	p(\varrho)  = a\varrho^\gamma, 
	\quad \Pi (\varrho) = a \frac{\varrho^\gamma}{\gamma -1} ,
	\quad a = \operatorname{const} > 0,
	\quad \gamma > 1 
\end{aligned}
$$

"""
