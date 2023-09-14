import streamlit as st
from PIL import Image
import os.path


dirpath = os.path.dirname(__file__)
for i in range(2):
    dirpath = os.path.dirname(dirpath)

imgs = os.path.join(dirpath, "figs/plots")

fig4_P1 = Image.open(os.path.join(imgs, "fig4_P1.png"))
fig4_P2 = Image.open(os.path.join(imgs, "fig4_P2.png"))
fig4_P3 = Image.open(os.path.join(imgs, "fig4_P3.png"))

fig8_P1 = Image.open(os.path.join(imgs, "fig8_P1.png"))
fig8_P2 = Image.open(os.path.join(imgs, "fig8_P2.png"))
fig8_P3 = Image.open(os.path.join(imgs, "fig8_P3.png"))

fig11_P1 = Image.open(os.path.join(imgs, "fig11_P1.png"))
fig11_P2 = Image.open(os.path.join(imgs, "fig11_P2.png"))
fig11_P3 = Image.open(os.path.join(imgs, "fig11_P3.png"))

r"""
# Лагранжевы элементы высших порядков

### Нелинейная задача (плотность на прямой $x_2 = 0$ в разные моменты времени)
"""

cols = st.columns(3)

with cols[0]:
    r"""
    ##### Лагранжевы элементы 1-й степени
    """
    st.image(fig4_P1)

with cols[1]:
    r"""
    ##### Лагранжевы элементы 2-й степени
    """
    st.image(fig4_P2)

with cols[2]:
    r"""
    ##### Лагранжевы элементы 3-й степени
    """
    st.image(fig4_P3)


r"""
Точечная линия — $\tau = 0,01$, пунктирная — $\tau = 0,005$, сплошная — $\tau = 0,0025$

Наиболее монотонное решение получается при элементах первой степени 

При $\tau = 0,01$ решение остается монотонным с увеличением степени конечных элементов

### Итерационная схема (плотность на прямой $x_2 = 0$ в разные моменты времени)
"""

cols = st.columns(3)

with cols[0]:
    r"""
    ##### Лагранжевы элементы 1-й степени
    """
    st.image(fig8_P1)

with cols[1]:
    r"""
    ##### Лагранжевы элементы 2-й степени
    """
    st.image(fig8_P2)

with cols[2]:
    r"""
    ##### Лагранжевы элементы 3-й степени
    """
    st.image(fig8_P3)

r"""
Пунктирная линия — $K = 1$, точечная — $K = 2$, сплошная — $K = 5$

При $K = 1$ решение получается немонотонным.

При остальных рассматриваемых $K$ решение остается монотонным с ростом степени конечных элементов.

### Динамика полной механической энергии
"""

cols = st.columns(3)

with cols[0]:
    r"""
    ##### Лагранжевы элементы 1-й степени
    """
    st.image(fig11_P1)

with cols[1]:
    r"""
    ##### Лагранжевы элементы 2-й степени
    """
    st.image(fig11_P2)

with cols[2]:
    r"""
    ##### Лагранжевы элементы 3-й степени
    """
    st.image(fig11_P3)

r"""
Пунктирная линия — $K = 1$, сплошная — $K = 5$

Уменьшение временного шага приводит к повышению точности выполнения закона сохранения

Использование конечных элементов 2-й степени приводит к уменьшению энергии итерационной схемы при $K=1$

При использовании конечных элементов 3-й степени и $\tau = 0.01$ итерационная схема неустойчива
"""