import streamlit as st


menu = st.sidebar.radio('***',
                        (
                            "Закон сохранения массы",
                            "Закон сохранения импульса",
                            "Закон сохранения полной механической энергии",
                        )
                        )

if menu == "Закон сохранения массы":
    r"""
    ## Закон сохранения массы
    
    Непосредственным интегрированием уравнения неразрывности по области $\Omega$ с учетом граничного условия получается
    закон сохранения массы
    
    $\begin{aligned}
    m (t) = m(0),\quad m(t) = \int_{\Omega}  \varrho(x, t) \ dx
    \end{aligned}$
    """

if menu == "Закон сохранения импульса":
    r"""
    ## Закон сохранения импульса

    Интегрируя уравнение движения по $\Omega$, получим
    
    $\begin{aligned}
    \int_{\Omega} \frac{\partial }{\partial t} (\varrho \bm u) \ dx + \int_{\partial \Omega} p(\varrho) \bm n \ dx = 0
    \end{aligned}$
    
    Тем самым имеет место
    
    $\begin{aligned}
    \bm I (t) = \bm I(0) - \int_{\partial \Omega} p(\varrho) \bm n \ dx, \quad \bm I(t) =  \int_{\Omega} \varrho \bm u \ dx
    \end{aligned}$
    """

if menu == "Закон сохранения полной механической энергии":
    r"""
    ## Закон сохранения полной механической энергии
    
    $\begin{aligned}
    E(t) = E(0), \quad E(t) =  \int_{\Omega}\left ( \frac{1}{2} \varrho |\bm u |^2 + \Pi(\varrho) \right ) \ dx
    \end{aligned}$
    
    Потенциал давления $\Pi(\varrho)$ определяется из уравнения
    
    $\begin{aligned}
    \varrho \frac{d \Pi}{d \varrho} - \Pi(\varrho) = p(\varrho)
    \end{aligned}$
    
    где
    
    $\begin{aligned}
    p(\varrho)  = a\varrho^\gamma, \quad \Pi (\varrho) = a \frac{\varrho^\gamma}{\gamma -1} ,\quad a = \operatorname{const} > 0,\quad \gamma > 1 
    \end{aligned}$
    """
