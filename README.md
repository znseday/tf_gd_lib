## tf_gd_lib [![Build Status](https://travis-ci.org/znseday/tf_gd_lib.svg?branch=master)](https://travis-ci.org/znseday/tf_gd_lib)

## Tabulated Function And Gradient Descent Math Lib

This library was written by Sergey Tsynikin as a coursework for OTUS courses.
This library is published under the Boost Software License, Version 1.0. So please make sure you agree with all terms of it before start using the library.
See accompanying file LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

##### In Russian language
Данная библиотека написана в качестве курсовой работы (курсового проекта) в рамках прохождения курса OTUS "Разработчик C++".
Автор библиотеки: Цыникин Сергей.
Библиотека опубликована под лицензией Boost Software License, Version 1.0. Пожалуйста, убедитесь, что Вы принимаете и согласны со всеми положениями данной лицензии перед тем, как начать использовать библиотеку и/или ее код тем или иным способом.
Текст лицензии см. в файле LICENSE_1_0.txt или по ссылке https://www.boost.org/LICENSE_1_0.txt

### Description

This small math library provides with solutions for two typical problems/tasks.

1. A flyweight functional class that works with tabulated functions.
2. A parameter optimization class that uses a gradient descent method. 

### Working with tabulated functions
The library provides with a flyweight class (TableFunction) that keeps a tabulated function, but allows to be treated as a typical continuous function.
This class uses linear interpolation/extrapolation and cubic spline calculations to get function value at any point.

Also, this class contains a few methods for common tasks, for example, to sort, to clear, to kill duplicates, to get min/max values, etc.
Anyway, internal data must be sorted before get started to use.

This class uses binary search for random access (a value at any point) that has logarithmic complexity O(log(n)).
However, the last point which be used are cached, so sequences access has constant time complexity O(1) in most of cases and linear complexity ( O(n-i) or O(i) ) in the worst cases where i is an index of a cached point.
A sequences access can be used only with linear interpolation and extrapolation. A sequences access can be used only with linear interpolation and extrapolation.

Cubic spline calculations can be used only for random access.

This class has operator(), and can be used as a callable object. In this case, only random access can be used.

### Parameter optimization using gradient descent method

The class GradDescent solves a problem of parameter optimization. This class works with TableFunction class for experimental data and with continuous one-variable function as a target function. However, amount of function parameters is unlimited.
This class is developed to work with any functions. It doesn't matter how sharp ravine of parameter surface is.
Such a good behavior is achieved by using independent descent rates that are altered automatically during the calculations.
For example, this class was tested for damped oscillations, linear, polynomial functions, Gaussian distributions, and any sums of these functions.
Also, this class contains a callback function for tracking a calculation process or stopping calculations at any time.

### Tests
The file tests.cpp contains typical examples of using the library.
