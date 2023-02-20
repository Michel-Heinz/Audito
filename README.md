# Audito
## General Information
Audito is an automatic differentiation tool for fortran. The mathematical operators have been overloaded to work with the newly defined types, which include not only the function value, but also the gradient, Hessian and Laplacian. This repository contains 3 modules, FG_m, FGL_m, and FGH_m, where FG calculates the value and gradient; FGL the value, gradient, and Laplacian; and FGH the value, gradient, and Hessian. You can use which ever you want, they do not depend on each other!
## How to use the code?
In order to use this code simply use the fortran module in your main code! Due to some compiler bugs, the code currently only works with the latest ifx and ifort compilers aswell as gfortran. If another compiler is used, a memory leak occurs, which can crash your PC!
## Examples
3 Examples are included, one for the FG module (value and gradient), one for the FGL module (value, gradient, and Laplacian), and one for the FGH module (value, gradient, and Hessian). The Example code is the optimization of Lennard-Jones clusters and it is commented. If you wish to use the modules for your own code, you will understand how to use it by reading the comments in the examples codes. In order to compile the examples use the following command (compiler: ifort, ifx, or gfortran; module: FG_m.f90, FGL_m.f90, FGH_m.f90; example: ExampleFG.f90, ExampleFGL.f90, ExampleFGH.f90; mathlibrary: -mkl, -llapack):
```
<compiler> <module> <example> <mathlibrary>
```
## Remarks
This program is work in progress! If there are any bugs or questions about the code and/or how to use it, feel free to write me!

