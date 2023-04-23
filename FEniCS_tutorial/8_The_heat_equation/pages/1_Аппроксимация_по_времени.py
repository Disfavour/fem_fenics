import streamlit as st


r"""
# Аппроксимация по времени

Номер временного слоя обозначим верхним индексом, а параметр дискретизации по времени - $\tau$;
тогда время на $n$-м временном слое:

$$
t^n = n \tau
$$

Таким образом, справедливол следующее:

$$
f(x, t^n) = f^n(x)
$$

Дискретизация с конечной разностью во времени сначала состоит из выбора краевой задачи на некотором временном слое:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

- \div \grad y^{n+1} + {\left( \frac {\partial y} {\partial t} \right)}^{n+1} = f^{n+1}
$$

Производная по времени может быть аппроксимирована (разностная производная назад):

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

{\left( \frac {\partial y} {\partial t} \right)}^{n+1} \approx \frac {y^{n+1} - y^{n}} {\tau}
$$

Получаем:

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

-\div \grad y^{n+1} + \frac {y^{n+1} - y^{n}} {\tau} \approx f^{n+1}
$$

Это дискретная по времени версия уравнения теплопроводности (неявный метод Эйлера).
"""
