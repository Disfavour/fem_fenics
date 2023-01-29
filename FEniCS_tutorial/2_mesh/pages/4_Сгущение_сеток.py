import streamlit as st
from fenics import *
import matplotlib.pyplot as plt


r"""
# Сгущение сеток

## refine(mesh, redistribute=true)

Создание более густой равномерной сетки.

### Параметры
- mesh (Mesh) - сетка.
- redistribute (bool) - перераспределять сетку, если она распределенная.

"""

with st.expander("Визуализация"):
    "### Параметры прямоугольной сетки `RectangleMesh`"

    col1, col2 = st.columns(2)

    with col1:
        "##### p1"
        p1 = Point(st.slider("x1", -10.0, 10.0, 0.0, key=11), st.slider("x2", -10.0, 10.0, 0.0, key=401))

    with col2:
        "##### p2"
        p2 = Point(st.slider("x1", -10.0, 10.0, 1.0, key=13), st.slider("x2", -10.0, 10.0, 1.0, key=402))

    col1, col2 = st.columns(2)

    with col1:
        nx1 = st.slider("nx1", 1, 100, 10, key=403)
        nx2 = st.slider("nx2", 1, 100, 10, key=404)
        diagonal = st.selectbox("diagonal", ("left", "right", "right/left", "left/right", "crossed"))

    with col2:
        number_refines = st.slider("Количество применений `refine`", 0, 10, 1)

    mesh = RectangleMesh(p1, p2, nx1, nx2, diagonal)

    col1, col2 = st.columns(2)

    with col1:
        plot(mesh)
        st.pyplot(plt.gcf())
        plt.clf()

        f"Количество ячеек ${mesh.num_cells()}$"

    with col2:
        for i in range(number_refines):
            mesh = refine(mesh)

        plot(mesh)
        st.pyplot(plt.gcf())
        plt.clf()

        f"Количество ячеек ${mesh.num_cells()}$"
