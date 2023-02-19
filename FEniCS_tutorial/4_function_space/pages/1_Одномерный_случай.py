import streamlit as st
import imageio
from os.path import abspath, join, dirname
from fenics import *
import matplotlib.pyplot as plt
import numpy as np

img_dir = join(dirname(dirname(abspath(__file__))), "images")

numerated_interval = imageio.imread_v2(join(img_dir, "numerated_interval.png"))
numerated_interval_3p = imageio.imread_v2(join(img_dir, "numerated_interval_3p.png"))

r"""
# Одномерный случай

Область:

$$
\Omega = \{x \ | \ l_0 \le x \le l_1\}.
$$

Узлы сетки:

$$
x_i = \frac {l_1 - l_0} {n} i, \quad i = 0, 1, \dots, n.
$$

Ячейки сетки:

$$
\Omega_i = \{ x \  | \  x_{i-1} \le x \le x_i  \}, \quad i = 1, 2, \dots, n.
$$

Сетка при $l_0 = 0, \ l_1 = 1$ и $n = 5$:
"""

plot(UnitIntervalMesh(5))
st.pyplot(plt.gcf())
plt.clf()

r"""
Конечный элемент представляет собой ячейку вместе с аппроксимацией на ней.

## Кусочно-линейная аппроксимация

Рассмотрим ячейку $[x_0, x_1]$.
"""

st.image(numerated_interval)

r"""
Аппроксимация по 2-м узлам полиномом 1-й степени

$$
f(x) = ax + b
$$

получается из системы:

$$
\begin{cases}
f_0 = a x_0 + b
\\
f_1 = a x_1 + b
\end{cases}.
$$

Искомый полином:

$$
\frac {f_0 - f_1} {x_0 - x_1} x + \frac {f_1 x_0 - f_0 x_1} {x_0 - x_1}.
$$

### Конечно-элементный базис

Линейно независимые (пробные) функции

$$
\varphi_i(x_j) = 
\begin{cases}
1, \  i = j
\\
0, \  i \ne j
\end{cases}
$$

равны единице в одном определённом узле и обращаются в ноль во всех других узлах.

Базисные функции непрерывны на ячейке.

Условия для базисных функций:

$$
\begin{cases}
\varphi_0(x_0) = 1
\\
\varphi_0(x_1) = 0
\end{cases},

\quad \quad

\begin{cases}
\varphi_1(x_0) = 0
\\
\varphi_1(x_1) = 1
\end{cases}.
$$

Базисные функции могут быть найдены из систем:

$$
\begin{cases}
a_0 x_0 + b_0 = 1
\\
a_0 x_1 + b_0 = 0
\end{cases},

\quad \quad

\begin{cases}
a_1 x_0 + b_1 = 0
\\
a_1 x_1 + b_1 = 1
\end{cases}.
$$

Базисные функции:

$$
\varphi_0(x) = \frac {1} {x_0 - x_1} x + \frac {x_1} {x_1 - x_0},
\\[0.3 cm]
\varphi_1(x) = \frac {1} {x_1 - x_0} x + \frac {x_0} {x_0 - x_1}.
$$

При $x_0 = 0$ и $x_1 = 1$ базисные функции принимают вид:
"""

col = st.columns(2)

with col[0]:
    r"""
    $$
    \varphi_0(x) = -x + 1,
    \\[0.3 cm]
    \varphi_1(x) = x.
    $$
    """
with col[1]:
    n = 10
    x = np.linspace(-1, 2, n)

    y_0 = -x + 1
    y_1 = x

    plt.plot(x, y_0)
    plt.plot(x, y_1)
    plt.grid()
    plt.legend([r"$\varphi_0$", r"$\varphi_1$"])
    st.pyplot(plt.gcf())
    plt.clf()

r"""
#### Глобальный базис
"""

with st.columns(2)[0]:
    n = 300
    x = np.linspace(-1, 2, n)

    y_0 = np.array(x)
    y_0[:n // 3] = x[n // 3: 2 * n // 3]
    y_0[n // 3: 2 * n // 3] = -x[n // 3: 2 * n // 3] + 1
    y_0[2 * n // 3:] = np.zeros(n // 3)

    y_1 = np.array(x)
    y_1[:n // 3] = np.zeros(n // 3)
    y_1[n // 3: 2 * n // 3] = x[n // 3: 2 * n // 3]
    y_1[2 * n // 3:] = -x[n // 3: 2 * n // 3] + 1

    plt.plot(x, y_0)
    plt.plot(x, y_1)
    plt.plot((0, 0), (-0.2, 1.2), "r--")
    plt.plot((1, 1), (-0.2, 1.2), "r--")
    plt.grid()
    plt.legend([r"$\varphi_0$", r"$\varphi_1$"])
    st.pyplot(plt.gcf())
    plt.clf()

r"""
## Кусочно-квадратичная аппроксимация

Рассмотрим ячейку $[x_0, x_2], \  x_0 < x_1 < x_2$.
"""

st.image(numerated_interval_3p)

r"""
Аппроксимация по 3-м узлам полиномом 2-й степени

$$
f(x) = a x^2 + b x + c.
$$

Применяем интерполяционный многочлен Лагранжа:

$$
L(x) = \sum \limits_{i=0}^n y_i l_i(x),
\\[0.3 cm]
l_i(x) = \prod \limits_{j=0, \ j \ne i}^n \frac {x - x_j} {x_i - x_j}.
$$

Искомый полином:

$$
f_0 \frac {(x - x_1)(x - x_2)} {(x_0 - x_1)(x_0 - x_2)}
+ f_1 \frac {(x - x_0)(x - x_2)} {(x_1 - x_0)(x_1 - x_2)}
+ f_2 \frac {(x - x_0)(x - x_1)} {(x_2 - x_0)(x_2 - x_1)}.
$$

### Конечно-элементный базис

Условия для базисных функций:

$$
\begin{cases}
\varphi_0(x_0) = 1
\\
\varphi_0(x_1) = 0
\\
\varphi_0(x_2) = 0
\end{cases},

\quad \quad

\begin{cases}
\varphi_1(x_0) = 0
\\
\varphi_1(x_1) = 1
\\
\varphi_1(x_2) = 0
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

С помощью интерполяционного многочлена Лагранжа получаем базисные функции:

$$
\varphi_0 = \frac {(x - x_1)(x - x_2)} {(x_0 - x_1)(x_0 - x_2)},
\\[0.3 cm]
\varphi_1 = \frac {(x - x_0)(x - x_2)} {(x_1 - x_0)(x_1 - x_2)},
\\[0.3 cm]
\varphi_2 = \frac {(x - x_0)(x - x_1)} {(x_2 - x_0)(x_2 - x_1)}.
$$

При $x_0 = 0, \ x_1 = 1, \ x_2 = 2$ базисные функции принимают вид:
"""

col = st.columns(2)

with col[0]:
    r"""
    $$
    \varphi_0(x) = \frac {1} {2} x^2 - \frac {3} {2} x + 1,
    \\[0.3 cm]
    \varphi_1(x) = - x^2 + 2 x,
    \\[0.3 cm]
    \varphi_2(x) = \frac {1} {2} x^2 - \frac {1} {2} x.
    $$
    """
with col[1]:
    n = 200
    x = np.linspace(0, 2, n)

    y_0 = 0.5 * x ** 2 - 1.5 * x + 1
    y_1 = - x ** 2 + 2 * x
    y_2 = 0.5 * x ** 2 - 0.5 * x

    plt.plot(x, y_0)
    plt.plot(x, y_1)
    plt.plot(x, y_2)
    plt.grid()
    plt.legend([r"$\varphi_0$", r"$\varphi_1$", r"$\varphi_2$"])
    st.pyplot(plt.gcf())
    plt.clf()

r"""
#### Глобальный базис
"""

with st.columns(2)[0]:
    n = 300
    x = np.linspace(-2, 4, n)

    y_0 = np.array(x)
    y_0[:n//3] =0.5 * x[n//2:5*n//6] ** 2 - 1.5 * x[n//2:5*n//6] + 1
    y_0[n//3:2*n//3] = 0.5 * x[n//3:2*n//3] ** 2 - 1.5 * x[n//3:2*n//3] + 1
    y_0[2*n//3:] = 0

    y_1 = np.array(x)
    y_1[:n // 3] = 0
    y_1[n // 3: 2 * n // 3] = - x[n // 3:2 * n // 3] ** 2 + 2 * x[n // 3: 2 * n // 3]
    y_1[2 * n // 3:] = 0

    y_2 = np.array(x)
    y_2[:n // 3] = 0
    y_2[n // 3:2 * n // 3] = 0.5 * x[n // 3:2 * n // 3] ** 2 - 0.5 * x[n // 3:2 * n // 3]
    y_2[2 * n // 3:] = 0.5 * x[n // 6: n // 2] ** 2 - 0.5 * x[n // 6: n // 2]

    plt.plot(x, y_0)
    plt.plot(x, y_1)
    plt.plot(x, y_2)
    plt.plot((0, 0), (-0.2, 1.2), "r--")
    plt.plot((2, 2), (-0.2, 1.2), "r--")
    plt.grid()
    plt.legend([r"$\varphi_0$", r"$\varphi_1$", r"$\varphi_2$"])
    st.pyplot(plt.gcf())
    plt.clf()
