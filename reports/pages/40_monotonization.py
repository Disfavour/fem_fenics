import streamlit as st


r'''
## Монотонизация уравнения переноса
$$
\frac {\partial u} {\partial t} + \frac {\partial u} {\partial x} = 0
$$


### Искусственная вязкость
$$
b = |\operatorname{grad} u|^\gamma = \lvert \frac {\partial u} {\partial x} \rvert ^\gamma
$$

$$
\int \limits_\Omega \frac {u^{n+1} - u^{n}} \tau u_t \ dx
+ \int \limits_\Omega \frac {\partial u^{n+\frac 1 2}} {\partial x} u_t \ dx
+ \int \limits_\Omega a \tau^2 b \frac {\partial u} {\partial x} \frac {\partial u_t} {\partial x} \ dx
= 0
$$

### Новая схема с $\kappa$
$$
\int \limits_\Omega \frac {u^{n+1} - u^{n}} \tau u_t \ dx
+ \int \limits_\Omega \frac {\partial u^{n+\frac 1 2}} {\partial x} u_t \ dx
+ \int \limits_\Omega (\kappa - \frac 1 2) \tau \frac {\partial} {\partial x} \left( \frac {u^{n+1} - u^{n}} \tau \right) u_t \ dx
= 0
$$
'''
