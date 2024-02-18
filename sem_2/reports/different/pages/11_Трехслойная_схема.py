import streamlit as st
from PIL import Image
from os.path import dirname, join, exists

dir = join(dirname(dirname(__file__)), 'images', 'err_etc')
hl, hr = 10, 1

r'''
# Трехслойная схема

$$
\frac {\partial \varphi} {\partial t} = \theta \frac {\varphi^{n+1} - \varphi^n} {\tau} + (1 - \theta) \frac {\varphi^{n} - \varphi^{n-1}} {\tau}
\\[0.5 cm]
\varphi^{n+\sigma} = \varphi^{n+1}
\\[0.5 cm]
$$

$$
\begin{aligned}
	& \theta \frac {h^{n+1} - h^n} {\tau} + (1 - \theta) \frac {h^{n} - h^{n-1}} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} u^{n+\sigma}) = 0
	\\[0.5 cm]
	& \theta \frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} + (1 - \theta) \frac {h^{n} u^{n} - h^{n-1} u^{n-1}} {\tau} + \frac {\partial} {\partial x} (h^{n+\sigma} (u^{n+\sigma})^2) + g h^{n+\sigma} \frac {\partial h^{n+\sigma}} {\partial x} = 0
\end{aligned}
$$
'''
