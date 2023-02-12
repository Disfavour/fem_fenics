import streamlit as st
from fenics import *
import matplotlib.pyplot as plt
import imageio
from os.path import abspath, join, dirname
import numpy as np

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
"""

col = st.columns(2)

with col[0]:
    r"""
    $$
    \varphi_0(x, y) = 1 - x - y,
    \\[0.3 cm]
    \varphi_1(x, y) = x,
    \\[0.3 cm]
    \varphi_2(x, y) = y.
    $$
    """

x = np.linspace(0, 1, 20)
y = np.array(x)

x, y = np.meshgrid(x, y)

z_0 = 1 - x - y
z_1 = x
z_2 = y

z = np.array(x)

for i in range(z_0.shape[0]):
    for j in range(z_0.shape[1]):
        z[i][j] = max(z_0[i][j], z_1[i][j], z_2[i][j])

with col[1]:
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.plot_surface(x, y, z)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.view_init(60, -60)

    st.pyplot(fig)
    plt.clf()

col = st.columns(3)

with col[0]:
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="3d")

    ax.plot_surface(x, y, z_0)

    ax.set_xlabel("x", fontsize=20)
    ax.set_ylabel("y", fontsize=20)
    ax.set_zlabel("z", fontsize=20)

    ax.set_title(r"$\varphi_0$", fontsize=20)

    ax.view_init(15, 120)

    st.pyplot(fig)
    plt.clf()

with col[1]:
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="3d")

    ax.plot_surface(x, y, z_1)

    ax.set_xlabel("x", fontsize=20)
    ax.set_ylabel("y", fontsize=20)
    ax.set_zlabel("z", fontsize=20)

    ax.set_title(r"$\varphi_1$", fontsize=20)

    ax.view_init(30, -120)

    st.pyplot(fig)
    plt.clf()

    cur_x = np.linspace(0, 1, 5)
    cur_z_1 = cur_x

    plt.plot(cur_x, cur_z_1)

    plt.xlabel("$x$", fontsize=20)
    plt.ylabel("$z$", fontsize=20)

    plt.suptitle(r"$\varphi_1$", fontsize=20)

    st.pyplot(plt.gcf())
    plt.clf()

with col[2]:
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection="3d")

    ax.plot_surface(x, y, z_2)

    ax.set_xlabel("x", fontsize=20)
    ax.set_ylabel("y", fontsize=20)
    ax.set_zlabel("z", fontsize=20)

    ax.set_title(r"$\varphi_2$", fontsize=20)

    st.pyplot(fig)
    plt.clf()

    cur_y = np.linspace(0, 1, 5)
    cur_z_2 = cur_y

    plt.plot(cur_y, cur_z_2)

    plt.xlabel("$y$", fontsize=20)
    plt.ylabel("$z$", fontsize=20)

    plt.suptitle(r"$\varphi_2$", fontsize=20)

    st.pyplot(plt.gcf())
    plt.clf()
