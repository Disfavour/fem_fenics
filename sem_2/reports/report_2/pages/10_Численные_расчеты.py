import streamlit as st
from PIL import Image
from os.path import dirname, join

images = join(dirname(dirname(__file__)), 'images')

r'''
## Численные расчеты

#### Динамика плотности в центре и экстремальных значений
'''
cols = st.columns(2)
with cols[0]:
    t1 = st.select_slider(r'Шаг по времени $\tau$', (0.01, 0.005, 0.0025), key=1)
with cols[1]:
    m1 = st.select_slider(r'Размер сетки', (100, 200), key=2)

cols = st.columns(3)
with cols[0]:
    r'''
    $$
    \sigma=0.5
    $$
    '''
    image = Image.open(join(images, f'fig2_w_s0.5_tau{t1}_ms{m1}.png'))
    st.image(image)
with cols[1]:
    r'''
    $$
    \sigma=0.75
    $$
    '''
    image = Image.open(join(images, f'fig2_w_s0.75_tau{t1}_ms{m1}.png'))
    st.image(image)
with cols[2]:
    r'''
    $$
    \sigma=1
    $$
    '''
    image = Image.open(join(images, f'fig2_w_s1_tau{t1}_ms{m1}.png'))
    st.image(image)

r'''
#### Плотность на различные моменты времени
'''
m2 = st.select_slider(r'Размер сетки', (100, 200), key=3)

cols = st.columns(3)
with cols[0]:
    r'''
    $$
    \sigma=0.5
    $$
    '''
    image = Image.open(join(images, f'fig4_w_s0.5_ms{m2}.png'))
    st.image(image)
with cols[1]:
    r'''
    $$
    \sigma=0.75
    $$
    '''
    image = Image.open(join(images, f'fig4_w_s0.75_ms{m2}.png'))
    st.image(image)
with cols[2]:
    r'''
    $$
    \sigma=1
    $$
    '''
    image = Image.open(join(images, f'fig4_w_s1_ms{m2}.png'))
    st.image(image)

r'''
- Линия из точек соответствует $\tau = 0.01$
- Штриховая линия соответствует $\tau = 0.005$
- Сплошная линия соответствует $\tau = 0.0025$

#### Динамика полной механической энергии при различных шагах по времени
'''
m3 = st.select_slider(r'Размер сетки', (100, 200), key=4)

cols = st.columns(3)
with cols[0]:
    r'''
    $$
    \sigma=0.5
    $$
    '''
    image = Image.open(join(images, f'fig11_w_s0.5_ms{m3}.png'))
    st.image(image)
with cols[1]:
    r'''
    $$
    \sigma=0.75
    $$
    '''
    image = Image.open(join(images, f'fig11_w_s0.75_ms{m3}.png'))
    st.image(image)
with cols[2]:
    r'''
    $$
    \sigma=1
    $$
    '''
    image = Image.open(join(images, f'fig11_w_s1_ms{m3}.png'))
    st.image(image)

r'''
#### Динамика полной механической энергии при различных $\sigma$
'''
cols = st.columns(2)
with cols[0]:
    t4 = st.select_slider(r'Шаг по времени $\tau$', (0.01, 0.005, 0.0025), key=5)
with cols[1]:
    m4 = st.select_slider(r'Размер сетки', (100, 200), key=6)

image = Image.open(join(images, f'fig_dif_sigma_w_t{t4}_ms{m4}.png'))
st.image(image)
