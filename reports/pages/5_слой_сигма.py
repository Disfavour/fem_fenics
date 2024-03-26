import streamlit as st

r'''
## $\varphi^{n+\sigma}$

$$
\begin{aligned}
    &h^{n+\sigma} = \sigma h^{n+1} + (1-\sigma) h^n
    \\[0.5 cm]
    &h^{n+1} = \frac {h^{n+\sigma} - (1-\sigma) h^n} {\sigma}
    \\[0.5 cm]
    &\frac {h^{n+1} - h^n} {\tau} = \frac {h^{n+\sigma} - (1-\sigma) h^n - \sigma h^n} {\tau \sigma}
    = \frac {h^{n+\sigma} - h^n} {\tau \sigma}
    \\[0.5 cm]
    &\frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} = \frac 1 \tau \left( \frac {h^{n+\sigma} - (1-\sigma) h^n} {\sigma}
    \frac {u^{n+\sigma} - (1-\sigma) u^n} {\sigma} - h^n u^n \right)
    \\[0.5 cm]
    &= \frac 1 {\tau \sigma^2} (h^{n+\sigma} u^{n+\sigma} - h^{n+\sigma} (1-\sigma) u^n - (1-\sigma) h^n u^{n+\sigma} + (1-\sigma) h^n (1-\sigma) u^n - \sigma^2 h^n u^n)
    \\[0.5 cm]
    &= \frac {h^{n+\sigma} u^{n+\sigma} - (1-\sigma) h^{n+\sigma} u^n - (1-\sigma) h^n u^{n+\sigma} + (1-\sigma)^2 h^n u^n - \sigma^2 h^n u^n} {\tau \sigma^2}
    \\[0.5 cm]
    &= \frac {h^{n+\sigma} u^{n+\sigma} - (1-\sigma) h^{n+\sigma} u^n - (1-\sigma) h^n u^{n+\sigma} + (1-2\sigma) h^n u^n} {\tau \sigma^2}
    \\[0.5 cm]
    &\frac {h^{n+1} u^{n+1} - h^n u^n} {\tau} = \frac {h^{n+\sigma} u^{n+\sigma} - h^n u^n} {\tau \sigma}
\end{aligned}
$$
'''
