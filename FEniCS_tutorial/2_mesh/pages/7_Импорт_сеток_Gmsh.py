import streamlit as st


r"""
# Импорт сеток Gmsh

Созданную сетку необходимо экспортировать в формате `.msh` со следующими параметрами:

- Version 2 ASCII Format
- Save all elements

Далее необходимо конвертировать сетку в поддерживаемый формат (`dolfin-convert`) и импортировать из файла обычным образом.
"""
