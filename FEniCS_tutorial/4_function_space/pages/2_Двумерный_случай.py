import streamlit as st
from fenics import *
import matplotlib.pyplot as plt
import imageio
from os.path import abspath, join, dirname
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

img_dir = join(dirname(dirname(abspath(__file__))), "images")

numerated_triangle = imageio.imread_v2(join(img_dir, "numerated_triangle.png"))

r"""
# Двумерный случай

Область:

$$
\Omega = \{x, \ y \ | \ l_0 \le x \le l_1, \ l_2 \le y \le l_3 \}.
$$

Узлы сетки:

$$
\begin{cases}
\displaystyle x_i = \frac {l_1 - l_0} {n} i, \quad i = 0, 1, \dots, n.
\\[0.3 cm]
\displaystyle y_j = \frac {l_3 - l_2} {m} j, \quad j = 0, 1, \dots, m.
\end{cases}.
$$

Ячейки сетки:

$$
\Omega_i = \{ x, y \  | \  x_{i-1} \le x \le x_i, \ y_{j-1} \le y \le y_j, \ y > x \ or \ y < x  \}
$$

Сетка при $l_0 = l_2 = 0, \ l_1 = l_3 = 1$ и $n = m = 5$:
"""

with st.columns(2)[0]:
    plot(UnitSquareMesh(5, 5, "left"))
    st.pyplot(plt.gcf())
    plt.clf()

r"""
## Кусочно-линейная аппроксимация

Рассмотрим ячейку:
"""

with st.columns(3)[1]:
    st.image(numerated_triangle)

r"""
Аппроксимация по 3-м узлам

$$
f(x, y) = ax + b y + c
$$

получается из системы:

$$
\begin{cases}
f_0 = a x_0 + b y_0 + c
\\
f_1 = a x_1 + b y_1 + c
\\
f_2 = a x_2 + b y_2 + c
\end{cases}.
$$

Искомый полином:

$$
\frac {(f_1 - f_2)(y_0 - y_1) - (f_0 - f_1)(y_1 - y_2)} {(x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2)} x
\\[0.3 cm]
+ \frac {(f_0 - f_1)((x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2))
- (x_0 - x_1)((f_1 - f_2)(y_0 - y_1) - (f_0 - f_1)(y_1 - y_2))}
{(y_0 - y_1)((x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2))} y
\\[0.3 cm]
+ \frac {f_0 (y_0 - y_1)((x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2))
- ((f_1 - f_2)(y_0 - y_1) - (f_0 - f_1)(y_1 - y_2))(y_0 - y_1)x_0

- ((f_0 - f_1)((x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2))
- (x_0 - x_1)((f_1 - f_2)(y_0 - y_1) - (f_0 - f_1)(y_1 - y_2)))y_0}

{(y_0 - y_1)((x_1 - x_2)(y_0 - y_1) - (x_0 - x_1)(y_1 - y_2))}
$$

### Конечно-элементный базис

Условия для базисных функций:

$$
\begin{cases}
\varphi_0(x_0, y_0) = 1
\\
\varphi_0(x_1, y_1) = 0
\\
\varphi_0(x_2, y_2) = 0
\end{cases},

\quad \quad

\begin{cases}
\varphi_1(x_0, y_0) = 0
\\
\varphi_1(x_1, y_1) = 1
\\
\varphi_1(x_2, y_2) = 0
\end{cases},

\quad \quad

\begin{cases}
\varphi_2(x_0) = 0
\\
\varphi_2(x_1) = 0
\\
\varphi_2(x_2) = 1
\end{cases}.
$$

Базисные функции могут быть найдены из систем:

$$
\begin{cases}
a_0 x_0 + b_0 y_0 + c_0 = 1
\\
a_0 x_1 + b_0 y_1 + c_0 = 0
\\
a_0 x_2 + b_0 y_2 + c_0 = 0
\end{cases},

\quad \quad

\begin{cases}
a_1 x_0 + b_1 y_0 + c_1 = 0
\\
a_1 x_1 + b_1 y_1 + c_1 = 1
\\
a_1 x_2 + b_1 y_2 + c_1 = 0
\end{cases},

\quad \quad

\begin{cases}
a_2 x_0 + b_2 y_0 + c_2 = 0
\\
a_2 x_1 + b_2 y_1 + c_2 = 0
\\
a_2 x_2 + b_2 y_2 + c_2 = 1
\end{cases}.
$$

Базисные функции:

$$
\varphi_0(x, y) = \frac {y_1 - y_2} {(x_0 - x_1)(y_1 - y_2) - (x_1 - x_2)(y_0 - y_1)} x
\\[0.3 cm]
- \frac {x_1 - x_2} {(x_0 - x_1)(y_1 - y_2) - (x_1 - x_2)(y_0 - y_1)} y
\\[0.3 cm]
- \frac {(y_1 - y_2) x_1 + (x_1 - x_2) y_1} {(x_0 - x_1)(y_1 - y_2) - (x_1 - x_2)(y_0 - y_1)},
\\[0.9 cm]
\varphi_1(x, y) = \frac {y_0 - y_2} {(x_1 - x_0)(y_0 - y_2) - (x_0 - x_2)(y_1 - y_0)} x
\\[0.3 cm]
- \frac {x_0 - x_2} {(x_1 - x_0)(y_0 - y_2) - (x_0 - x_2)(y_1 - y_0)} y
\\[0.3 cm]
- \frac {(y_0 - y_2) x_0 + (x_0 - x_2) y_0} {(x_1 - x_0)(y_0 - y_2) - (x_0 - x_2)(y_1 - y_0)},
\\[0.9 cm]
\varphi_2(x, y) = \frac {y_1 - y_0} {(x_2 - x_1)(y_1 - y_0) - (x_1 - x_0)(y_2 - y_1)} x
\\[0.3 cm]
- \frac {x_1 - x_0} {(x_2 - x_1)(y_1 - y_0) - (x_1 - x_0)(y_2 - y_1)} y
\\[0.3 cm]
- \frac {(y_1 - y_0) x_1 + (x_1 - x_0) y_1} {(x_2 - x_1)(y_1 - y_0) - (x_1 - x_0)(y_2 - y_1)}.
$$

При 

$$
\begin{cases}
(x_0, y_0) = (0, 0)
\\
(x_1, y_1) = (1, 0)
\\
(x_2, y_2) = (0, 1)
\end{cases}
$$

базисные функции принимают вид:

$$
\varphi_0(x, y) = 1 - x - y,
\\[0.3 cm]
\varphi_1(x, y) = x,
\\[0.3 cm]
\varphi_2(x, y) = y.
$$
"""

col = st.columns(3)

with col[0]:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    v = np.array([
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
    ])

    ax.scatter3D(v[:, 0], v[:, 1], v[:, 2])

    verts = [
        [v[0], v[1], v[2]],
    ]

    ax.add_collection3d(
        Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))

    plt.xlabel("$x$")
    plt.ylabel("$y$")
    fig.suptitle(r"$\varphi_0$")

    #ax.view_init(30, -120)

    st.pyplot(fig)

with col[1]:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    v = np.array([
        [0, 0, 0],
        [0, 1, 0],
        [1, 0, 1],
    ])

    ax.scatter3D(v[:, 0], v[:, 1], v[:, 2])

    verts = [
        [v[0], v[1], v[2]],
    ]

    ax.add_collection3d(
        Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))

    plt.xlabel("$x$")
    plt.ylabel("$y$")
    fig.suptitle(r"$\varphi_1$")

    ax.view_init(30, -120)

    st.pyplot(fig)

with col[2]:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    v = np.array([
        [0, 0, 0],
        [0, 1, 1],
        [1, 0, 0],
    ])

    ax.scatter3D(v[:, 0], v[:, 1], v[:, 2])

    verts = [
        [v[0], v[1], v[2]],
    ]

    ax.add_collection3d(
        Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))

    plt.xlabel("$x$")
    plt.ylabel("$y$")
    fig.suptitle(r"$\varphi_2$")

    ax.view_init(30, -120)

    st.pyplot(fig)

"#### Глобальный базис"

with st.columns(2)[0]:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    v = np.array([
        [-1, 0, 0],
        [-1, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 1, 0],
        [1, -1, 0],
        [1, 0, 0],

        [-2, -2, 0],
        [-2, 2, 0],
        [2, 2, 0],
        [2, -2, 0],
    ])

    ax.scatter3D(v[:, 0], v[:, 1], v[:, 2])

    verts = [
        [v[0], v[1], v[3]],
        [v[1], v[4], v[3]],
        [v[0], v[2], v[3]],
        [v[4], v[6], v[3]],
        [v[2], v[5], v[3]],
        [v[3], v[5], v[6]],
        [v[7], v[8], v[9], v[10]]
    ]

    ax.add_collection3d(
        Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25))

    plt.xlabel("$x$")
    plt.ylabel("$y$")

    st.pyplot(fig)
