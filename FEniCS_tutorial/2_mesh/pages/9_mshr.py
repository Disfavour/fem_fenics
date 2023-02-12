import streamlit as st
from fenics import *
import mshr
import matplotlib.pyplot as plt


r"""
# mshr

mshr - это компонент FEniCS для генерации сеток.

Может быть установлен:
`sudo apt-get install python3-mshr`

Cначала создается область (domain), затем создается сетка (mesh).

## Создание сетки

mshr.generate_mesh(domain, resolution, backend="cgal")

- domain (CSGGeometry) - область.
- resolution (float) - разрешение сетки. Размер ячейки будет приблизительно равен диаметру домена, деленному на
  разрешение.
- backend (string) - используемый внутренний интерфейс.

## Создание областей

### Двумерный многоугольник

mshr.Polygon(points)

- points (List[Point]) - список точек, описывающих двумерный многоугольник.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=901)
        domain = mshr.Polygon([Point(0, 0), Point(1, 0), Point(0, 1)])
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Окружность

mshr.Circle(c, r, segments=32)

- c (Point) - центр.
- r (float) - радиус.
- segments (int) - количество сегментов при вычислении полигональной аппроксимации.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=902)
        c = Point(0, 0)
        r = st.slider("r", 0.01, 10.0, 1.0, key=903)
        segments = st.slider("segments", 0, 100, 5, key=904)
        domain = mshr.Circle(c, r, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Прямоугольник

mshr.Rectangle(p1, p2)

- p1 (Point) - вершина.
- p2 (Point) - вершина.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=905)
        domain = mshr.Rectangle(Point(0, 0), Point(2, 1))
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Эллипс

mshr.Ellipse(c, a, b, segments=32)

- c (Point) - центр.
- a (float) - горизонтальная полуось.
- b (float) - вертикальная полуось.
- segments (int) - количество сегментов при вычислении полигональной аппроксимации.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=906)
        c = Point(0, 0)
        a = st.slider("a", 0.01, 10.0, 2.0, key=907)
        b = st.slider("b", 0.01, 10.0, 1.0, key=908)
        segments = st.slider("segments", 0, 100, 5, key=909)
        domain = mshr.Ellipse(c, a, b, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Цилиндр

mshr.Cylinder(top, bottom, top_radius, bottom_radius, segments=32)

- top (Point) - центр в верхней части цилиндра.
- bottom (Point) - центр в нижней части цилиндра.
- top_radius (float) - радиус верхней части цилиндра.
- bottom_radius (float) - радиус в нижней части цилиндра.
- segments (int) - количество ограничивающих граней при генерации многогранного приближения.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=910)
        top = Point(0, 0, 1)
        bottom = Point(0, 0, 0)
        top_radius = st.slider("top_radius", 0.01, 10.0, 2.0, key=911)
        bottom_radius = st.slider("bottom_radius", 0.01, 10.0, 1.0, key=912)
        segments = st.slider("segments", 0, 100, 5, key=913)
        domain = mshr.Cylinder(top, bottom, top_radius, bottom_radius, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Прямоугольный параллелепипед

mshr.Cylinder(p1, p2)

- p1 (Point) - вершина.
- p2 (Point) - вершина.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=914)
        p1 = Point(0, 0, 0)
        p2 = Point(1, 2, 3)
        domain = mshr.Box(p1, p2)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Конус

mshr.Cone(top, bottom, r, segments=32)

- top (Point) - центр в верхней части конуса.
- bottom (Point) - центр в нижней части конуса.
- r (float) - верхний радиус конуса.
- segments (int) - количество ограничивающих граней при генерации многогранного приближения.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=915)
        top = Point(0, 0, 1)
        bottom = Point(0, 0, 0)
        r = st.slider("r", 0.01, 10.0, 2.0, key=916)
        segments = st.slider("segments", 0, 100, 5, key=917)
        domain = mshr.Cone(top, bottom, r, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Сфера

mshr.Sphere(center, radius, segments=10)

- center (Point) - центр сферы.
- radius (float) - радиус сферы.
- segments (int) - разрешение при генерации многогранного приближения.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 5.0, key=918)
        center = Point(0, 0, 0)
        radius = st.slider("radius", 0.01, 10.0, 2.0, key=919)
        segments = st.slider("segments", 0, 100, 5, key=920)
        domain = mshr.Sphere(center, radius, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Тетраэдр

mshr.Tetrahedron(a, b, c, d)

- a (Point) - вершина.
- b (Point) - вершина.
- c (Point) - вершина.
- d (Point) - вершина.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 10.0, 2.0, key=921)
        a = Point(0, 0, 0)
        b = Point(1, 0, 0)
        c = Point(0, 1, 0)
        d = Point(0, 0, 1)
        domain = mshr.Tetrahedron(a, b, c, d)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Эллипсоид

mshr.Ellipsoid(center, a, b, c, segments=15)

- center (Point) - центр.
- a (float) - полуось в направлении x.
- b (float) - полуось в направлении y.
- c (float) - полуось в направлении z.
- segments (int) - разрешение при генерации многогранного приближения.
"""

with st.expander("Пример"):
    with st.echo():
        resolution = st.slider("resolution", 0.01, 20.0, 7.0, key=922)
        center = Point(0, 0, 0)
        a = st.slider("a", 0.01, 10.0, 1.0, key=923)
        b = st.slider("b", 0.01, 10.0, 2.0, key=924)
        c = st.slider("c", 0.01, 10.0, 1.0, key=925)
        segments = st.slider("segments", 0, 100, 50, key=926)
        domain = mshr.Ellipsoid(center, a, b, c, segments)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Триангулярная 3-мерная поверхность из файла

mshr.Surface3D(filename, degenerate_tolerance=1e-12)

- filename (string) - имя файла (путь). Поддерживаются off, stl, vtp файлы.
- degenerate_tolerance (float) - допуск ухудшения.

## Операции над областями

### Объединение `+`
"""

with st.expander("Пример"):
    with st.echo():
        domain = mshr.Circle(Point(-0.8, 0), 1) + mshr.Circle(Point(0.8, 0), 1)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Пересечение `*`
"""

with st.expander("Пример"):
    with st.echo():
        domain = mshr.Circle(Point(-0.8, 0), 1) * mshr.Circle(Point(0.8, 0), 1)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()

r"""
### Разница `-`
"""

with st.expander("Пример"):
    with st.echo():
        domain = mshr.Circle(Point(-0.8, 0), 1) - mshr.Circle(Point(0.8, 0), 1)
        mesh = mshr.generate_mesh(domain, resolution)

        plot(mesh)

    with st.columns(2)[0]:
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        st.pyplot(plt.gcf())
        plt.clf()
