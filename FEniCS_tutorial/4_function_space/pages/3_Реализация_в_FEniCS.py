import streamlit as st
from os.path import abspath, join, dirname
import imageio


img_dir = join(dirname(dirname(abspath(__file__))), "images")

lagrange_interval = imageio.imread_v2(join(img_dir, "d1m.png"))
lagrange_triangle = imageio.imread_v2(join(img_dir, "lagrange2d.png"))
lagrange_tetrahedron = imageio.imread_v2(join(img_dir, "lagrange3d.png"))


r"""
# Функциональные пространства в FEniCS

## FunctionSpace(mesh, family, degree, form_degree=None, constrained_domain=None, restriction=None)

Создание конечно-элементного функционального пространства.

### Параметры
- mesh (Mesh) - сетка.
- family (string) - тип конечных элементов.
  - "CG" - Lagrange
  - "DG" - Discontinuous Lagrange
  - "ARG" - Argyris
  - "AW" - Arnold-Winther
  - "BDFM" - Brezzi-Douglas-Fortin-Marini
  - "BDM" - Brezzi-Douglas-Marini
  - "B" - Bubble
  - "CR" - Crouzeix-Raviart
  - "DRT" - Discontinuous Raviart-Thomas
  - "HER" - Hermite
  - "MTW" - Mardal-Tai-Winther
  - "MOR" - Morley
  - "N1curl" - Nedelec 1st kind H(curl)
  - "N2curl" - Nedelec 2nd kind H(curl)
  - "Q" - Quadrature
  - "RT" - Raviart-Thomas
  - "R" - Real
- degree (int) - степень элемента.
- form_degree (int) - степень формы (обозначение FEEC, используемое, когда поле рассматривается как k-форма).
- constrained_domain - ограниченный поддомен с функцией отображения.
- restriction - ограничение элемента (например, на грани ячейки).
"""

with st.expander("Пример использования"):
    """
    Дискретное функциональное пространство над единичным квадратом:
    ```python
    mesh = UnitSquare(32,32)
    V = FunctionSpace(mesh, "CG", 1)
    ```
    """

r"""
## 1-мерные лагранжевы элементы

Для более высокой степени элемента задействуются точки внутри отрезка.
"""

st.image(lagrange_interval)

r"""
## 2-мерные лагранжевы элементы

Для более высокой степени элемента задействуются точки на ребрах и внутри треугольника.
"""
st.image(lagrange_triangle)

r"""
## 3-мерные лагранжевы элементы

Для более высокой степени элемента задействуются точки на ребрах и внутри тетраэдра.
"""
st.image(lagrange_tetrahedron)
