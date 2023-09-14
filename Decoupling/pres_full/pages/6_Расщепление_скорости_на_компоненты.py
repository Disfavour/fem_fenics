import streamlit as st


r"""
# Расщепление скорости на компоненты

Пусть $\bm u = (u, v)$, где $u$ - компонента по $x_1$, $v$ - компонента по $x_2$

### Расщепление типа Якоби

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\frac{\varrho_{n+1}^{k+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1}^{k+1} \bm u_{n+1}^k)
    = 0
    \\[0.5 cm]
    &\frac{\varrho_{n+1}^{k+1} u_{n+1}^{k+1} - \varrho_{n} u_{n}} {\tau}
    + \frac {\partial (u_{n+1}^{k} \varrho_{n+1}^{k+1} u_{n+1}^{k+1})} {\partial x_1}
    + \frac {\partial (v_{n+1}^{k} \varrho_{n+1}^{k+1} u_{n+1}^{k+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1}^{k+1})} {\partial x_1}
    = 0
    \\[0.5 cm]
    &\frac{\varrho_{n+1}^{k+1} v_{n+1}^{k+1} - \varrho_{n} v_{n}} {\tau}
    + \frac {\partial (u_{n+1}^{k} \varrho_{n+1}^{k+1} v_{n+1}^{k+1})} {\partial x_1}
    + \frac {\partial (v_{n+1}^{k} \varrho_{n+1}^{k+1} v_{n+1}^{k+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1}^{k+1})} {\partial x_2}
    = 0
\end{aligned}
$$

### Расщепление типа Зейделя

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\frac{\varrho_{n+1}^{k+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1}^{k+1} \bm u_{n+1}^k)
    = 0
    \\[0.5 cm]
    &\frac{\varrho_{n+1}^{k+1} u_{n+1}^{k+1} - \varrho_{n} u_{n}} {\tau}
    + \frac {\partial (u_{n+1}^{k} \varrho_{n+1}^{k+1} u_{n+1}^{k+1})} {\partial x_1}
    + \frac {\partial (v_{n+1}^{k} \varrho_{n+1}^{k+1} u_{n+1}^{k+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1}^{k+1})} {\partial x_1}
    = 0
    \\[0.5 cm]
    &\frac{\varrho_{n+1}^{k+1} v_{n+1}^{k+1} - \varrho_{n} v_{n}} {\tau}
    + \frac {\partial (\textcolor{red} {u_{n+1}^{k+1}} \varrho_{n+1}^{k+1} v_{n+1}^{k+1})} {\partial x_1}
    + \frac {\partial (v_{n+1}^{k} \varrho_{n+1}^{k+1} v_{n+1}^{k+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1}^{k+1})} {\partial x_2}
    = 0
\end{aligned}
$$
"""

with st.expander("Линеаризованная схема"):
    r"""
### Линеаризованная схема с разделением скорости на компоненты

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\frac {\varrho_{n+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1} \bm u_{n})
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} u_{n+1} - \varrho_{n} u_{n}} {\tau}
    + \div(\bm u_{n} \varrho_{n+1} u_{n+1})
    + \frac {\partial p(\varrho_{n+1})} {\partial x_1}
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} v_{n+1} - \varrho_{n} v_{n}} {\tau}
    + \div(\bm u_{n} \varrho_{n+1} v_{n+1})
    + \frac {\partial p(\varrho_{n+1})} {\partial x_2}
    = 0
\end{aligned}
$$


##### Расщепление типа Якоби

Преобразуем, чтобы явно выделить компоненты с прошлого временного слоя

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\frac {\varrho_{n+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1} \bm u_{n})
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} u_{n+1} - \varrho_{n} u_{n}} {\tau}
    + \frac {\partial (u_n \varrho_{n+1} u_{n+1})} {\partial x_1}
    + \frac {\partial (v_n \varrho_{n+1} u_{n+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1})} {\partial x_1}
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} v_{n+1} - \varrho_{n} v_{n}} {\tau}
    + \frac {\partial (u_n \varrho_{n+1} v_{n+1})} {\partial x_1}
    + \frac {\partial (v_n \varrho_{n+1} v_{n+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1})} {\partial x_2}
    = 0
\end{aligned}
$$

##### Расщепление типа Зейделя

На основе предыдущего преобразования можно предложить при расчете $v_{n+1}$ использовать полученное из предыдущего
уравнения $u_{n+1}$ вместо $u_n$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
    &\frac {\varrho_{n+1} - \varrho_{n}} {\tau}
    + \div(\varrho_{n+1} \bm u_{n})
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} u_{n+1} - \varrho_{n} u_{n}} {\tau}
    + \frac {\partial (u_n \varrho_{n+1} u_{n+1})} {\partial x_1}
    + \frac {\partial (v_n \varrho_{n+1} u_{n+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1})} {\partial x_1}
    = 0
    \\[0.5 cm]
    &\frac {\varrho_{n+1} v_{n+1} - \varrho_{n} v_{n}} {\tau}
    + \frac {\partial (u_{n+1} \varrho_{n+1} v_{n+1})} {\partial x_1}
    + \frac {\partial (v_n \varrho_{n+1} v_{n+1})} {\partial x_2}
    + \frac {\partial p(\varrho_{n+1})} {\partial x_2}
    = 0
\end{aligned}
$$
"""
