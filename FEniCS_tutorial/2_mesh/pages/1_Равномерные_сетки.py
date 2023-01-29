import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Равномерные сетки

Для создания некоторых сеток используется объект представляющий точку - Point.

### Point(x1=0.0, x2=0.0, x3=0.0)

Создание точки $(x_1, x_2, x_3)$.
"""

r"""
### IntervalMesh(n, a, b)

Создание равномерной сетки на отрезке $[a, b]$.

#### Параметры

- n (int) - количество ячеек.

- a (float) - левая граница отрезка.

- b (float) - правая граница отрезка.
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        n = st.slider("n", 1, 100, 10)
        a, b = st.slider("a, b", -2.0, 2.0, (0.0, 1.0))

    with col2:
        plot(IntervalMesh(n, a, b))
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### UnitIntervalMesh(n)

Аналог IntervalMesh, но только для отрезка $[0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitIntervalMesh(n)
```

```python
IntervalMesh(n, 0, 1)
```
"""

r"""
### RectangleMesh(p1, p2, nx1, nx2, diagonal="right")

Создание треугольной равномерной сетки на прямоугольнике.

#### Параметры

- p1 (Point) - точка с минимальными значениями $x_1, x_2$ прямоугольника.

- p2 (Point) - точка с максимальными значениями $x_1, x_2$ прямоугольника.

- nx1 (int) - количество ячеек по $x_1$.

- nx2 (int) - количество ячеек по $x_2$.

- diagonal (str) - строка, указывающая направление диагоналей ячеек.

  Может принимать следующие значения:

    - "left"
    - "right"
    - "right/left"
    - "left/right"
    - "crossed"
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        "##### p1"
        p1 = Point(st.slider("x1", -10.0, 10.0, 0.0, key=11), st.slider("x2", -10.0, 10.0, 0.0, key=12))

    with col2:
        "##### p2"
        p2 = Point(st.slider("x1", -10.0, 10.0, 1.0, key=13), st.slider("x2", -10.0, 10.0, 1.0, key=14))

    col1, col2 = st.columns(2)

    with col1:
        nx1 = st.slider("nx1", 1, 100, 10, key=15)
        nx2 = st.slider("nx2", 1, 100, 10, key=16)
        diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"))

    with col2:
        plot(RectangleMesh(p1, p2, nx1, nx2, diagonal))
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### UnitSquareMesh(x1n, x2n, diagonal="right")

Аналог RectangleMesh, но только для квадрата $[0, 1] \times [0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitSquareMesh(x1n, x2n, diagonal)
```

```python
RectangleMesh(Point(0, 0), Point(1, 1), x1n, x2n, diagonal)
```
"""

r"""
### BoxMesh(p1, p2, nx1, nx2, nx3)

Создание равномерной сетки из тетраэдров на прямоугольном параллелепипеде.

#### Параметры

- p1 (Point) - точка с минимальными значениями $x_1, x_2, x_3$ прямоугольного параллелепипеда.

- p2 (Point) - точка с максимальными значениями $x_1, x_2, x_3$ прямоугольного параллелепипеда.

- nx1 (int) - количество ячеек по $x_1$.

- nx2 (int) - количество ячеек по $x_2$.

- nx3 (int) - количество ячеек по $x_3$.
"""

with st.expander("Визуализация"):
    col1, col2 = st.columns(2)

    with col1:
        "##### p1"
        p1 = Point(
            st.slider("x1", -10.0, 10.0, 0.0, key=21),
            st.slider("x2", -10.0, 10.0, 0.0, key=22),
            st.slider("x3", -10.0, 10.0, 0.0, key=23),
        )

    with col2:
        "##### p2"
        p2 = Point(
            st.slider("x1", -10.0, 10.0, 1.0, key=24),
            st.slider("x2", -10.0, 10.0, 1.0, key=25),
            st.slider("x3", -10.0, 10.0, 1.0, key=26),
        )

    col1, col2 = st.columns(2)

    with col1:
        nx1 = st.slider("nx1", 1, 100, 10, key=27)
        nx2 = st.slider("nx2", 1, 100, 10, key=28)
        nx3 = st.slider("nx3", 1, 100, 10, key=29)

    with col2:
        plot(BoxMesh(p1, p2, nx1, nx2, nx3))
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### UnitCubeMesh(nx1, nx2, nx3)

Аналог BoxMesh, но только для куба $[0, 1] \times [0, 1] \times [0, 1]$.

Следующие два вызова эквивалентны.

```python
UnitCubeMesh(nx1, nx2, nx3)
```

```python
BoxMesh(Point(0, 0, 0), Point(1, 1, 1), nx1, nx2, nx3)
```
"""
