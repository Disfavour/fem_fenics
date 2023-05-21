import streamlit as st


r"""
# Некоторые операции

## Скалярное произведение

Для векторов $a = (a_1, \dots, a_N)$, $b = (b_1, \dots, b_N)$:

$$
a \cdot b = \sum_{i=1}^N a_i b_i
$$

Для тензоров $A = \{A_{i,j}\}_{i,j=1}^N$, $B = \{B_{i,j}\}_{i,j=1}^N$:

$$
A : B = \sum_{i,j=1}^N A_{i,j} B_{i,j}
$$

## Тензорное произведение

Для векторов $a,\ b$:

$$
a \otimes b = \{ a_i b_j \}_{i,j=1}^N
$$

Пример тензорного произведения векторов:

$$
a \otimes b
=
\begin{bmatrix}a_1 \\ a_2 \\ a_3 \\ a_4\end{bmatrix}  
\begin{bmatrix}b_1 & b_2 & b_3\end{bmatrix} = 
\begin{bmatrix}
a_1b_1 & a_1b_2 & a_1b_3 \\
a_2b_1 & a_2b_2 & a_2b_3 \\
a_3b_1 & a_3b_2 & a_3b_3 \\
a_4b_1 & a_4b_2 & a_4b_3
\end{bmatrix}
$$

Пример тензорного произведения матриц:

$$
\begin{bmatrix}
    a_{1,1} & a_{1,2} \\
    a_{2,1} & a_{2,2} \\
  \end{bmatrix}
  \otimes
  \begin{bmatrix}
    b_{1,1} & b_{1,2} \\
    b_{2,1} & b_{2,2} \\
  \end{bmatrix}
  =
  \begin{bmatrix}
    a_{1,1} \begin{bmatrix}
      b_{1,1} & b_{1,2} \\
      b_{2,1} & b_{2,2} \\
    \end{bmatrix} & a_{1,2} \begin{bmatrix}
      b_{1,1} & b_{1,2} \\
      b_{2,1} & b_{2,2} \\
    \end{bmatrix} \\[0.5 cm]
    a_{2,1} \begin{bmatrix}
      b_{1,1} & b_{1,2} \\
      b_{2,1} & b_{2,2} \\
    \end{bmatrix} & a_{2,2} \begin{bmatrix}
      b_{1,1} & b_{1,2} \\
      b_{2,1} & b_{2,2} \\
    \end{bmatrix} \\
  \end{bmatrix}
  =
  \begin{bmatrix}
    a_{1,1} b_{1,1} & a_{1,1} b_{1,2} & a_{1,2} b_{1,1} & a_{1,2} b_{1,2} \\
    a_{1,1} b_{2,1} & a_{1,1} b_{2,2} & a_{1,2} b_{2,1} & a_{1,2} b_{2,2} \\
    a_{2,1} b_{1,1} & a_{2,1} b_{1,2} & a_{2,2} b_{1,1} & a_{2,2} b_{1,2} \\
    a_{2,1} b_{2,1} & a_{2,1} b_{2,2} & a_{2,2} b_{2,1} & a_{2,2} b_{2,2} \\
  \end{bmatrix}
$$

## Векторное произведение

Для векторов $a = [a_1, a_2, a_3]$, $b = [b_1, b_2, b_3]$:

$$
a \times b = 
\begin{vmatrix}
i &j &k \\
a_1 &a_2 &a_3 \\
b_1 &b_2 &b_3
\end{vmatrix}
= (a_2 b_3 - a_3 b_2; \  a_3 b_1 - a_1 b_3; \  a_1 b_2 - a_2 b_1)
$$

## Внешнее произведение

Для векторов в индексной записи:

$$
(\mathbf{u} \otimes \mathbf{v})_{ij} = u_i v_j
$$

Внешнее произведение $u \otimes v$ эквивалентно матричному умножению $u v^T$.

Пример:

$$
\mathbf{u} \otimes \mathbf{v} = \mathbf{u}\mathbf{v}^\textsf{T} =
  \begin{bmatrix}u_1 \\ u_2 \\ u_3 \\ u_4\end{bmatrix}
    \begin{bmatrix}v_1 & v_2 & v_3\end{bmatrix} =
  \begin{bmatrix}
    u_1 v_1 & u_1 v_2 & u_1 v_3 \\
    u_2 v_1 & u_2 v_2 & u_2 v_3 \\
    u_3 v_1 & u_3 v_2 & u_3 v_3 \\
    u_4 v_1 & u_4 v_2 & u_4 v_3
  \end{bmatrix}.
$$
"""
