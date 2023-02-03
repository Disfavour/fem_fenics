import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Трансформация сеток

Совершая действия над координатами узлов сетки можно изменять её.

Для повышения точности вблизи определенной границы можно использовать уплотнение сетки в этом направлении.

Пусть есть сетка с равномерно распределенными узлами $x_0, x_1, \dots, x_n$ на отрезке [a, b].

Трансформация координат

$$
\begin{equation}
\xi = \frac {x - a} {b - a}
\end{equation}
$$

преобразует в область $[0, 1]$.

Отображение

$$
\begin{equation}
\eta = \xi^s
\end{equation}
$$

для некоторого $s > 1$ растягивает сетку в направлении $\xi = 0$ $(x = a)$, в то время как для $s < 1$ растяжение
происходит в направлении $\xi = 1$ $(x = b)$.
"""

with st.expander("Визуализация"):
    col = st.columns(2)
    with col[0]:
        s = st.slider("s", 0.1, 5.0, 1.5)

    n = 10

    mesh = UnitSquareMesh(n, n)

    x = mesh.coordinates()[:, 0]

    mesh.coordinates()[:, 0] = x ** s

    plot(mesh)

    with col[1]:
        st.pyplot(plt.gcf())
        plt.clf()

r"""
Отображение $\kappa$ от $\eta \in [0, 1]$ обратно к $[a, b]$ дает новые, растянутые координаты:

$$
\begin{equation}
\kappa = a + (b - a) \left( \frac {x - a} {b - a} \right)^s
\end{equation}
$$

Одним из способов создания более сложных геометрий является преобразование координат вершин в прямоугольной сетке в
соответствии с некоторой формулой.

Для создания фигуры, ограниченной двумя радиусами можно использовать стандартное отображение из полярных координат в
декартовы.

Есть прямоугольник $a \le x \le b$ и $0 \le y \le 1$, тогда отображение

$$
\begin{equation}
\tilde{x} = x \cos (\Theta y)
\end{equation}
$$

$$
\begin{equation}
\tilde{y} = x \sin (\Theta y)
\end{equation}
$$

трансформирует прямоугольник в часть круга $\Theta$ радиан.
"""

with st.expander("Реализация и визуализация"):
    with st.echo():
        from fenics import *
        import numpy as np
        import matplotlib.pyplot as plt

        Theta = pi / 2
        a, b = 1, 5.0
        nr = 20
        nt = 40
        mesh = RectangleMesh(Point(a, 0), Point(b, 1), nr, nt, "right")

        x = mesh.coordinates()[:, 0]
        y = mesh.coordinates()[:, 1]
        s = 1.5

        def denser(x, y):
            return [a + (b - a) * ((x - a) / (b - a)) ** s, y]

        x_bar, y_bar = denser(x, y)
        xy_bar_coor = np.array([x_bar, y_bar]).transpose()
        mesh.coordinates()[:] = xy_bar_coor

        def cylinder(r, s):
            return [r * np.cos(Theta * s), r * np.sin(Theta * s)]

        x_hat, y_hat = cylinder(x_bar, y_bar)
        xy_hat_coor = np.array([x_hat, y_hat]).transpose()
        mesh.coordinates()[:] = xy_hat_coor

        plot(mesh)

    col = st.columns(2)

    with col[0]:
        st.pyplot(plt.gcf())
        plt.clf()
