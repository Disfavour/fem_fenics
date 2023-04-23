import streamlit as st


r"""
# FEniCS реализация

Самый простой способ реализовать переменный коэффициент $\kappa(x)$ - это определить `Expression`, которое зависит от
координат.

```python
class K(Expression):
    def set_k_values(self, k_0, k_1):
        self.k_0, self.k_1 = k_0, k_1
        
    def eval(self, value, x):
        "Установка в value[0] значения в x"
        tol = 1E-14
        if x[1] <= 0.5 + tol:
            value[0] = self.k_0
        else:
            value[0] = self.k_1
            
            
# Инициализация
kappa = K()
kappa.set_k_values(1, 0.01)
```

Метод eval обеспечивает большую гибкость в определении функций, но недостатком является то, что FEniCS будет вызывать
eval в Python для каждого узла x, что является медленным.

Альтернативный метод заключается в использовании строкового выражения C++. Это можно сделать с помощью тернарного
оператора:

```python
tol = 1E-14
k_0 = 1.0
k_1 = 0.01
kappa = Expression("x[1] <= 0.5 + tol ? k_0 : k_1", degree=0, tol=tol, k_0=k_0, k_1=k_1)
```

Этот метод определения переменных коэффициентов работает, если подобласти представляют
собой простые формы, которые могут быть выражены в терминах геометрических неравенств.

Для более сложных подобластей нужно использовать маркирование ("SubDomain" и "MeshFunction").

Рассмотрим следующее определение границы x = 0:

```python
def boundary(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol)
```

Это определение границы на самом деле является сокращением к более общему SubDomain. SubDomain - это класс, который
определяет подобласть в пространстве в терминах метода, который возвращает True для точек, принадлежащих подобласти, и
False для точек, которые не принадлежат подобласти.

Вот как указать границу x = 0 в качестве подобласти:

```python
class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-14
        return on_boundary and near(x[0], 0, tol)
        
        
boundary = Boundary()
bc = DirichletBC(V, Constant(0), boundary)
```

Мы будем использовать два подкласса SubDomain для определения двух подобластей $\Omega_0$ и $\Omega_1$:

```python
tol = 1E-14

class Omega_0(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= 0.5 + tol
        
class Omega_1(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= 0.5 - tol
```

Вершины при $x_2=0.5$ должны принадлежать обеим подобластям, чтобы FEniCS корректно определял какой подобласти относится
каждая ячейка.

Чтобы определить переменный коэффициент $\kappa$, будем использовать MeshFunction. MeshFunction - дискретная функция,
которая может быть вычислена по набору так называемых объектов сетки. Объект сетки в FEniCS - это либо вершина, ребро,
грань, либо ячейка (треугольник или тетраэдр).
MeshFunction над ячейками подходит для представления подобластей,
в то время как MeshFunction над гранями используется для представления фрагментов внешних или внутренних границ.
MeshFunction над ячейками также может использоваться для представления граничных маркеров для уточнения сетки
(refine).

Для определения подобластей $\Omega_0$ и $\Omega_1$ используем MeshFunction. Первым аргументом конструктора является тип
маркера ("int", "size_t", "double", "bool"). В качестве аргументов можно указать имя файла для инициализации
MeshFunction.

Создадим MeshFunction с неотрицательными целочисленными значениями ("size_t"):

```python
materials = MeshFunction("size_t", mesh, mesh.topology().dim())
```

Далее мы используем две подобласти, чтобы промаркировать соответствующие ячейки:

```python
subdomain_0 = Omega_0()
subdomain_1 = Omega_1()
subdomain_0.mark(materials, 0)
subdomain_1.mark(materials, 1)
```

Это установит значения MeshFunction materials равными 0 для каждой ячейки, принадлежащей $\Omega_0$, и 1 для всех ячеек,
принадлежащих $\Omega_1$. В качестве альтернативы, мы можем использовать
следующий эквивалентный код для маркировки ячеек:

```python
materials.set_all(0)
subdomain_1.mark(materials, 1)
```

Чтобы проверить значения функции сетки и убедиться, что мы действительно
правильно определили наши поддомены, мы можем просто построить график функции сетки:

```python
plot(materials)
plt.show()
```

Возможно, мы также захотим сохранить значения MeshFunction для последующего использования:

```python
File("materials.xml") << materials
```

Получить MeshFunction из файла можно аналогичным образом:

```python
File("materials.xml") >> materials
```

Для использовать значения MeshFunction materials для определения
переменного коэффициента $\kappa$, мы создаем Expression:

```python
class K(UserExpression):
    def __init__(self, materials, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.materials = materials
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.materials[cell.index] == 0:
            values[0] = self.k_0
        else:
            values[0] = self.k_1
            

kappa = K(materials, k_0, k_1, degree=0)
```

Эта версия функции "eval" имеет дополнительный аргумент "cell", который мы можем использовать для соответствующих
проверок.

Поскольку мы используем геометрические тесты для определения двух подобластей для $\Omega_0$
и $\Omega_1$, метод MeshFunction может показаться ненужным усложнением простого метода, использующего выражение
с if-тестом. Однако, определение подобластей может быть доступно в виде MeshFunction
(из файла), сгенерированной как часть процесса генерации сетки,
а не как простой геометрический тест.

## Использование фрагментов кода C++ для определения подобластей

Классы SubDomain и Expression очень удобны, но
их использование приводит к вызовам функций из C++ в Python для каждого узла в
сетке. Поскольку это сопряжено со значительными затратами, необходимо использовать код на C++
, если возникает проблема с производительностью.

Вместо того, чтобы писать подкласс SubDomain на Python, мы можем вместо этого использовать
инструмент CompiledSubDomain в FEniCS, чтобы указать подобласть в
коде C++ и тем самым ускорить наш код. Рассмотрим определение классов
Omega_0 и Omega_1 выше на Python. Ключевые строки, определяющие эти
поддомены, могут быть выражены в синтаксисе C++ и переданы в качестве аргументов для
CompiledSubDomain следующим образом:

```python
tol = 1E-14
subdomain_0 = CompiledSubDomain("x[1] <= 0.5 + tol", tol=tol)
subdomain_1 = CompiledSubDomain("x[1] >= 0.5 - tol", tol=tol)
```

Результирующие объекты, subdomain_0 и subdomain_1, могут использоваться как обычные объекты SubDomain.

Скомпилированные строки поддомена также могут быть применены для указания границ:

```python
boundary_R = CompiledSubDomain("on_boundary && near(x[0], 1, tol)", tol=1E-14)
```

Также возможно передать строку C++ (без параметров) непосредственно
в качестве третьего аргумента DirichletBC без явного построения
CompiledSubDomain:

```python
bc1 = DirichletBC(V, value, "on_boundary && near(x[0], 1, tol)")
```
"""

with st.expander("Реализация"):
    r"""
    ```python
    from fenics import *
    import matplotlib.pyplot as plt
    
    tol = 1E-14
    k_0 = 1.0
    k_1 = 0.01

    mesh = UnitSquareMesh(20, 20)

    V = FunctionSpace(mesh, "P", 1)
    
    materials = MeshFunction("size_t", mesh, mesh.topology().dim())
    
    subdomain_0 = CompiledSubDomain("x[1] <= 0.5 + tol", tol=tol)
    subdomain_1 = CompiledSubDomain("x[1] >= 0.5 - tol", tol=tol)
    
    materials.set_all(0)
    subdomain_1.mark(materials, 1)
    
    class K(UserExpression):
        def __init__(self, materials, k_0, k_1, **kwargs):
            super().__init__(**kwargs)
            self.materials = materials
            self.k_0 = k_0
            self.k_1 = k_1
    
        def eval_cell(self, values, x, cell):
            if self.materials[cell.index] == 0:
                values[0] = self.k_0
            else:
                values[0] = self.k_1


    kappa = K(materials, k_0, k_1, degree=0)
        
    bcs = DirichletBC(V, Constant(0), "on_boundary")
    
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = kappa*dot(grad(u), grad(v))*dx
    L = f*v*dx
    
    u = Function(V)
    solve(a == L, u, bcs)
        
    plt.colorbar(plot(u))
    plot(mesh)
    plt.show()
    ```
    """
