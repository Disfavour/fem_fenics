import streamlit as st


r"""
## Консервативные схемы второго порядка аппроксимации по времени для задач гидродинамики в приближении мелкой воды

#### Проблема

Построение численных методов для решения задач
гидродинамики на современных компьютерах параллельной
архитектуры, которые обеспечивают устойчивость
приближенного решения и наследования основных свойств
дифференциальной задачи (положительность решения,
консервативность).

#### Актуальность

- Модели сжимаемой жидкости широко применяются
  в различных областях, что делает моделирование
  критически важным.
- Современные вычислительные ресурсы и методы
  позволяют более глубоко исследовать процессы,
  происходящие в сжимаемых жидкостях, и
  разрабатывать более точные и эффективные
  численные методы.
  
#### Новизна

- Использование новых переменных в задачах
  гидродинамики, которые обеспечивают
  неотрицательность плотности.
- Выполнение законов сохранения на дискретном
  уровне.
- Применение безусловно устойчивых (в линейном
  случае) схем второго порядка по времени.
- Удобная вычислительная реализация за счет
  развязывания отдельных уравнений системы.
  
#### Цель работы

Разработка вычислительного алгоритма второго порядка
аппроксимации по времени, его программная реализация и
верификация для численного решения задач гидродинамики
в приближении мелкой воды с наследованием
положительности плотности и законов сохранения на
дискретном уровне.

#### Задачи по диссертации

1. Математические модели мелкой воды
2. Аппроксимации по времени и пространству
3. Программная реализация для задач с плоским дном и
   без учета силы Кориолиса и демонстрация
   работоспособности
4. Численные эксперименты верификации
   вычислительного алгоритма и прикладного ПО на
   тестовых задачах
5. Реализация и результаты моделирования при
   использовании более общих моделей мелкой воды
   (рельеф дна, сила Кориолиса)
"""
