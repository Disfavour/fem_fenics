import streamlit as st


r"""
# Вариационная формулировка

Для использования FEniCS необходимо привести дифференциальное уравнение к вариационной формулировке.

### Приведение к вариационной формулировке

Умножаем уравнение Пуассона на **тестовую функцию** $v$ и интегрируем по $\Omega$:

$$
\begin{equation} \tag{4}
- \int \limits_\Omega (\operatorname{div} \operatorname{grad} u) v d \bm x = \int \limits_\Omega f v d \bm x.
\end{equation}
$$

Интегрируем по частям, чтобы избавиться от вторых производных от **пробной функциии** $u$:

$$
\begin{equation} \tag{5}
- \int \limits_\Omega (\operatorname{div} \operatorname{grad} u) v d \bm x = \int \limits_\Omega \operatorname{grad} (u)
\cdot \operatorname{grad} (v) d \bm x - 
\int \limits_{\partial\Omega} \frac {\partial u} {\partial \vec{n}} v \operatorname{ds},
\end{equation}
$$

где $\frac {\partial u} {\partial \vec{n}} = \nabla u \cdot \vec{n}$ является производной от $u$ во внешнем нормальном
направлении $\vec{n}$ на границе.

Тестовая функция $v$ обращается в ноль на тех участках границы, где известно решение $u$. В рамках рассматрвиаемой 
задачи это означает, что $v = 0$ на $\partial \Omega$:

$$
\begin{equation} \tag{6}
\int \limits_\Omega \operatorname{grad} (u) \cdot \operatorname{grad} (v) d \bm x = \int \limits_\Omega f v d \bm x.
\end{equation}
$$

Если потребовать, чтобы уравнение выполнялось для всех тестовых функций $v$ в **тестовом функциональном пространстве**
$\tilde{V}$, то получится четко определенная математическая задача, которая однозначно определяет решение $u$,
которое лежит в **пробном функциональном пространстве** $V$.

### Каноническая форма вариационной задачи

Найти $u \in V$ такое, что

$$
\begin{equation} \tag{7}
a(u, v) = l(v) \quad \forall v \in \tilde{V}.
\end{equation}
$$

Для уравнения Пуассона:

$$
\begin{equation} \tag{8}
a(u, v) = \int \limits_\Omega \operatorname{grad} (u) \cdot \operatorname{grad} (v) d \bm x \\
\end{equation}
$$

$$
\begin{equation} \tag{9}
l(v) = \int \limits_\Omega f v d \bm x.
\end{equation}
$$
"""
