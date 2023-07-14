import streamlit as st
from PIL import Image
import os.path


dirpath = os.path.dirname(__file__)
for i in range(2):
    dirpath = os.path.dirname(dirpath)

imgs = os.path.join(dirpath, "figs/plots")

nl_jac = Image.open(os.path.join(imgs, "nl_jac_P1.png"))
nl_sei = Image.open(os.path.join(imgs, "nl_sei_P1.png"))

r"""
# Сравнение расщеплений
"""

cols = st.columns(2)

with cols[0]:
    r"""
    ### Расщепление типа Якоби
    """
    st.image(nl_jac)

with cols[1]:
    r"""
    ### Расщепление типа Зейдель
    """
    st.image(nl_sei)

r"""
Относительно левой шкалы построена синяя кривая, которая соответствует энергии от времени нелинейной задачи

Относительно правой шкалы построены кривые, которые соответствуют модулю разности
указанного расщепления при $K$ итерациях и нелинейной задачи

Можно сделать вывод, что расщепление типа Зейдель не даёт весомых преимуществ, а расщепление Якоби может быть обратно
преобразовано к векторной форме, что дает возможность использовать векторные элементы
"""
