import streamlit as st


r"""
# Аппроксимаци по времени

Пусть $\tau$ - шаг равномерной, сетки во времени, такой, что $\varphi_n = \varphi(t_n), \ t_n = n \tau$,
$n = 0,1, ..., N, \ N\tau = T$.

При построении и исследовании аппроксимации по времени основное внимание
уделяется выполнению соответствующих законов сохранения (априорных оценок).

Для приближенного решения задачи используется чисто неявная 
схема. Решение на новом слое определяется в этом случае из 

$$
\begin{aligned}
 \frac{\varrho_{n+1} - \varrho_{n} }{\tau} + A(\bm u_{n+1}) \varrho_{n+1}  = 0,
\end{aligned}
$$

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}  }{\tau } + 
 A(\bm u_{n+1}) (\varrho_{n+1} \bm u_{n+1}) + \grad p(\varrho_{n+1}) = 0, 
 \quad n = 0, 1, ..., N-1 ,
\end{aligned} 
$$

Основные свойства приближенного решения связаны с выполнением законов сохранения
массы и полной энергии.
Будем считать, что для плотность положительна.

На каждом слое по времени $\varrho_{n}>0$, $n = 0, 1, ..., N$.

Интегрирование по области уравнения дает

$$
\begin{aligned}
 (\varrho_{n+1}, 1) = (\varrho_{n}, 1),
 \quad n = 0, 1, ..., N-1 .
\end{aligned} 
$$

Равенство - дискретный аналог закона сохранения массы (5).
Закону сохранения импульса (6) сопоставляется равенство

$$
\begin{aligned}
 (\varrho_{n+1} \bm u_{n+1}, 1) = (\varrho_{n} \bm u_{n}, 1)
 - \tau \int_{\partial \Omega} p(\varrho_{n+1})  \bm n d \bm x ,
 \quad n = 0, 1, ..., N-1 , 
\end{aligned} 
$$

Оценка для полной механической энергии устанавливается следующим образом.
Домножая уравнение на $\bm u_{n+1}$ и интегрируя по $\Omega$, получим

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\begin{split}
 \left ( \frac{\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}  }{\tau } , \bm u_{n+1} \right )  & + 
 (A(\bm u_{n+1}) (\varrho_{n+1} \bm u_{n+1}), \bm u_{n+1} ) \\
 & + (\grad p(\varrho_{n+1}), \bm u_{n+1} ) = 0 . 
\end{split}
\end{aligned}
$$

Для первого слагаемого имеем

$$
\begin{split}
 \frac{\varrho_{n+1} \bm u_{n+1} - \varrho_{n} \bm u_{n}  }{\tau } \bm u_{n+1} 
 & = \frac{1}{2} \frac{\varrho_{n+1} |\bm u_{n+1}|^2 - \varrho_{n} |\bm u_{n}|^2  }{\tau } \\
 & + \frac{1}{2} \frac{\varrho_{n+1} |\bm u_{n+1}|^2 + \varrho_{n} |\bm u_{n}|^2 - 2 \varrho_{n} \bm u_{n} \bm u_{n+1} }{\tau } \\
 & =  \frac{1}{2} \frac{\varrho_{n+1} |\bm u_{n+1}|^2 - \varrho_{n} |\bm u_{n}|^2  }{\tau } \\
 & + \frac{1}{2} \frac{\varrho_{n+1} - \varrho_{n} }{\tau }|\bm u_{n+1}|^2
 + \frac{1}{2} \varrho_{n}  \frac{ |\bm u_{n+1} - \bm u_{n}|^2 }{\tau } \\
 & \leq  \frac{1}{2} \frac{\varrho_{n+1} |\bm u_{n+1}|^2 - \varrho_{n} |\bm u_{n}|^2  }{\tau } 
 + \frac{1}{2} \frac{\varrho_{n+1} - \varrho_{n} }{\tau }|\bm u_{n+1}|^2 .
\end{split}
$$

В силу этого из уравнения следует

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
\begin{split}
 \frac{1}{2} \left ( \frac{\varrho_{n+1} |\bm u_{n+1}|^2 - \varrho_{n} |\bm u_{n}|^2  }{\tau } , 1 \right ) 
 & + \frac{1}{2} \left (\frac{\varrho_{n+1} - \varrho_{n} }{\tau }, |\bm u_{n+1}|^2 \right ) \\
 & + (A(\bm u_{n+1}) (\varrho_{n+1} \bm u_{n+1}), \bm u_{n+1} ) \\
 & + (\grad p(\varrho_{n+1}), \bm u_{n+1} ) \leq  0 . 
\end{split}
\end{aligned}
$$

Принимая во внимание (\ref{22}) и определение оператора $A$, имеем 

$$
\begin{split}
\frac{1}{2} \left (\frac{\varrho_{n+1} - \varrho_{n} }{\tau }, |\bm u_{n+1}|^2 \right ) 
 & + (A(\bm u_{n+1}) (\varrho_{n+1} \bm u_{n+1}), \bm u_{n+1} ) \\
 & = (A(\bm u_{n+1}) (\varrho_{n+1} \bm u_{n+1}), \bm u_{n+1} )
 - \frac{1}{2} (A(\bm u_{n+1}) \varrho_{n+1}), |\bm u_{n+1}|^2) = 0 .
\end{split} 
$$

Это дает возможность записать неравенство в виде

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \frac{1}{2} \left ( \frac{\varrho_{n+1} |\bm u_{n+1}|^2 - \varrho_{n} |\bm u_{n}|^2  }{\tau } , 1 \right ) 
 - (p(\varrho_{n+1}), \div \bm u_{n+1} ) \leq  0 . 
\end{aligned} 
$$

Для оценки второго слагаемого в левой части неравенства привлекается 
дискретный аналог (10).
Домножим уравнение неразрывности на ${\displaystyle \frac{d \varPi}{d \varrho}  (\varrho_{n+1})}$:

$$
\begin{aligned}
 \frac{\varrho_{n+1} - \varrho_{n} }{\tau}  \frac{d \varPi}{d \varrho} (\varrho_{n+1}) 
 + A(\bm u_{n+1}) \varrho_{n+1}  \frac{d \varPi }{d \varrho} (\varrho_{n+1}) = 0 . 
\end{aligned} 
$$

Имеет место равенство

$$
 \varPi(\varrho_{n+1}) - \varPi(\varrho_{n}) = \frac{d \varPi}{d \varrho} (\varrho_{n+1}) (\varrho_{n+1}-\varrho_{n})
 - \frac{1}{2} \frac{d^2 \varPi}{d \varrho^2} (\widetilde{\varrho}^{n+1}) (\varrho_{n+1}) (\varrho_{n+1}-\varrho_{n})^2 ,
$$

где

$$
 \widetilde{\varrho}^{n+1} \in [ \min(\varrho_{n}, \varrho_{n+1}), \ \max(\varrho_{n}, \varrho_{n+1}) ] .
$$

Пусть

$$
\begin{aligned}
 \frac{d^2 \varPi}{d \varrho^2} (\varrho) \geq 0 , 
\end{aligned} 
$$

При этих естественных ограничениях получим

$$
 \frac{d \varPi}{d \varrho} (\varrho_{n+1}) (\varrho_{n+1}-\varrho_{n}) \geq  
 \varPi(\varrho_{n+1}) - \varPi(\varrho_{n}) .
$$

Для второго слагаемого с учетом (8) имеем

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{split}
 A(\bm u_{n+1}) \varrho_{n+1}  \frac{d \varPi }{d \varrho} (\varrho_{n+1}) 
 & =  \varrho_{n+1}  \frac{d \varPi }{d \varrho} (\varrho_{n+1}) \div \bm u_{n+1} + \bm u_{n+1} \grad \varPi (\varrho ^{n+1}) \\
 & = \div(\varPi(\varrho_{n+1}) \bm u_{n+1}) + p(\varrho_{n+1}) \div \bm u_{n+1} .
\end{split} 
$$

С учетом этого интегрирование уравнения дает

$$
\def \grad {\operatorname{grad}}
\def \div {\operatorname{div}}
\def \rot {\operatorname{rot}}

\begin{aligned}
 \left (\frac{\varPi(\varrho_{n+1}) - \varPi(\varrho_{n})}{\tau} , 1 \right ) + (p(\varrho_{n+1}), \div \bm u_{n+1}) \leq 0 .  
\end{aligned} 
$$

Объединяя (\ref{28}) и (\ref{31}) придем к неравенству

$$
\begin{aligned}
 \left (\frac{1}{2} \varrho_{n+1} |\bm u_{n+1}|^2 + \varPi(\varrho_{n+1}) , 1 \right ) 
 \leq \left (\frac{1}{2} \varrho_{n} |\bm u_{n}|^2 + \varPi(\varrho_{n}) , 1 \right ) .
\end{aligned} 
$$

Сопоставляя уравнение с (12), можем сделать вывод о том, что на дискретном уровне вместо сохранения
полной энергии имеет место убывание энергии.

Для чисто неявной схемы,
которая дает приближенное решение задачи, 
выполнен закон сохранения массы, закон сохранения импульса, оценка для полной механической энергии.
"""
