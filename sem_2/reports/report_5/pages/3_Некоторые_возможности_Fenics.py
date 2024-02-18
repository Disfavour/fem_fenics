import streamlit as st
from fenics import *


'## Некоторые возможности Fenics'

'Координаты узлов сетки'
with st.echo():
    mesh = UnitIntervalMesh(4)
    mesh.coordinates()
st.code(mesh.coordinates())

'Значения функции в узлах сетки'
with st.echo():
    P = FunctionSpace(mesh, 'P', 1)
    u = Function(P)
    u.assign(project(Expression('x[0]', degree=1), P))
    u.compute_vertex_values()
st.code(u.compute_vertex_values())

'Значения функции в узлах интерполяции (degrees of freedom)'
st.code('u.vector()[:]')
st.code(u.vector()[:])
r'''
- `u.vector()[:]` - слайсинг по вектору, можно менять значения
- `u.vector().get_local()` - получения локальной копии (np.array)
- `u.vector().set_local(arg)` - установить значения arg (np.array) в вектор
'''

'Карта преобразования порядка узлов интерполяции в порядок узлов сетки'
st.code('dof_to_vertex_map(P)')
st.code(dof_to_vertex_map(P))

'Преобразование порядка узлов интерполяции в порядок узлов сетки'
st.code('u.vector()[dof_to_vertex_map(P)]')
st.code(u.vector()[dof_to_vertex_map(P)])

'Карта преобразования порядка узлов сетки в порядок узлов интерполяции'
st.code('vertex_to_dof_map(P)')
st.code(vertex_to_dof_map(P))

'Преобразование порядка узлов сетки в порядок узлов интерполяции'
st.code('u.compute_vertex_values()[vertex_to_dof_map(P)]')
st.code(u.compute_vertex_values()[vertex_to_dof_map(P)])

'### Сложные (векторные, mixed) функциональные пространства'
'Значения функции в узлах сетки'
with st.echo():
    H = FiniteElement('P', mesh.ufl_cell(), 1)
    U = FiniteElement('P', mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, MixedElement([H, U]))
    w = Function(W)
    w.assign(project(Expression(('x[0]', '2*x[0]'), degree=1), W))
    w.compute_vertex_values()
st.code(w.compute_vertex_values())
'Половины массива соответствуют подпространствам'

'Значения функции в узлах интерполяции'
st.code('w.vector()[:]')
st.code(w.vector()[:])

'Карта преобразования порядка узлов интерполяции в порядок узлов сетки'
st.code('dof_to_vertex_map(W)')
st.code(dof_to_vertex_map(W))

st.error('Неверное преобразование порядка степеней свободы в порядок узлов сетки')
st.code('w.vector()[dof_to_vertex_map(W)]')
st.code(w.vector()[dof_to_vertex_map(W)])

'Порядок узлов интерполяции подпространств'
st.code('W.sub(0).dofmap().dofs(), W.sub(1).dofmap().dofs()')
st.code((W.sub(0).dofmap().dofs(), W.sub(1).dofmap().dofs()))

'Значения функции в узлах интерполяции подпространств'
st.code('w.vector()[W.sub(0).dofmap().dofs()], w.vector()[W.sub(1).dofmap().dofs()]')
st.code((w.vector()[W.sub(0).dofmap().dofs()], w.vector()[W.sub(1).dofmap().dofs()]))

'Карта преобразования порядка узлов интерполяции в порядок узлов сетки для подпространств'
st.code('dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()], dof_to_vertex_map(W)[W.sub(1).dofmap().dofs()]')
st.code((dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()], dof_to_vertex_map(W)[W.sub(1).dofmap().dofs()]))

'Преобразование порядка узлов интерполяции в порядок узлов сетки для подпространств'
st.code('w.vector()[dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()]], w.vector()[dof_to_vertex_map(W)[W.sub(1).dofmap().dofs()]]')
st.code((w.vector()[dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()]], w.vector()[dof_to_vertex_map(W)[W.sub(1).dofmap().dofs()]]))

'### Конечно-элементные и сеточные функции'

'Cеточные функции на сетке функционального пространства'
st.code('mesh.coordinates().flatten()')
st.code(mesh.coordinates().flatten())
with st.echo():
    q1 = 10 + mesh.coordinates().flatten()
st.code(q1)
with st.echo():
    q2 = 1 + 3 * mesh.coordinates().flatten()
st.code(q2)

'Изменение функции путем замены значений вектора степеней свободы'
with st.echo():
    q = Function(W)
    q.vector()[dof_to_vertex_map(W)[W.sub(0).dofmap().dofs()]] = q1
    q.vector()[dof_to_vertex_map(W)[W.sub(1).dofmap().dofs()]] = q2
    q.compute_vertex_values()
st.code(q.compute_vertex_values())

'''
## Замечания

Это верно только для лагранжевых элементов 1-й степени

Для элементов более высокой степени количество степеней свободы превышает количество узлов и преобразование будет более сложным
'''
