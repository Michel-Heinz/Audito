# Audito
## General Information
Audito is an automatic differentiation tool for fortran. The mathematical operators have been overloaded to work with the newly defined types, which include not only the function value, but also the gradient, Hessian and Laplacian. This repository contains 3 modules, FG_m, FGL_m, and FGH_m, where FG calculates the value and gradient; FGL the value, gradient, and Laplacian; and FGH the value, gradient, and Hessian. You can use which ever you want, they do not depend on each other!
## What operations are currently implemented?
Let's define $x$ and $y$ as the new type FG_t or FGX_t (where X is either H or L) and $a$ to be any real (32 bit, 64 bit, or 128 bit) or any integer (8 bit, 16 bit, 32 bit, or 64 bit). The following operations are implemented:

| $x+y$ | $\text{e}^x$ | $a\cdot x$ |
| :---: | :---: | :---: |
| $x-y$ | $\sqrt{x}$ | $x\cdot a$ |
| $x\cdot y$ | $\text{ln}(x)$ |  $a + x$ |
|  $\frac{x}{y}$ | $\text{cos}(x)$ |  $x + a$ |
| $x^y$ | $\text{sin}(x)$ | $a - x$ |
|  | $\text{tan}(x)$ | $x - a$ |
|  | $-x$ | $\frac{a}{x}$ |
|  | $\text{abs}(x)$ | $\frac{x}{a}$ |
|  |  | $x^a$ |
|  |  | $a^x$ |


The FG_t or FGX_t types currently can only be defined for 64 bit reals! This will however be expanded in future releases.
## How does it work?
Overloading the operators enables a very easy implementation of functions! The comlicated part is the initialization of the variables. Consider the following example:

We have two atoms, each with $x$, $y$, and $z$ coordinates. I define a position vector $\vec{x}$ as follows:

$\vec{x} = (x_1, y_1, z_1, x_2, y_2, z_2)$

Let's define a function $f(\vec{x})$, which depends on the positions of the atoms. The Gradient of this function $\nabla f(\vec{x})$ looks like this:

$\nabla f(\vec{x}) = \left(\frac{f(\vec{x})}{x_1}, \frac{f(\vec{x})}{y_1}, \frac{f(\vec{x})}{z_1}, \frac{f(\vec{x})}{x_2}, \frac{f(\vec{x})}{y_2}, \frac{f(\vec{x})}{z_2}\right)^\text{T}$

Let's initialize the coordinates as an FG_t type. Each coordinate will have its function value x%f and its gradient x%g. Thus, $x_1$, $y_1$, and $y_2$ will look like this:

$x_1$%f = 0

$x_1$%g = $\left(1, 0, 0, 0, 0, 0\right)^\text{T}$

$y_1$%f = 0

$y_1$%g = $\left(0, 1, 0, 0, 0, 0\right)^\text{T}$

$y_2$%f = 0

$y_2$%g = $\left(0, 0, 0, 0, 1, 0\right)^\text{T}$

For each variable, its component in the gradient is set to 1 and the rest to zero. This is the result of the function $f(\vec{x})$, of which we only that it depends on these 6 variable (therefore the gradient of each variable has 6 components). Thus, we need to tell audito, the component of the gradient that belongs to each variable and set it to 1. The following fortran code initializes the gradients of each of the variables.
```
iniGrad = 0._r8 !Has size 6

do i = 1, SIZE(x)
    iniGrad(i) = 1._r8
    x(i) = FG_t(0._r8, iniGrad) !x has size 6 and represents the positons vector (type(FG_t) :: x(6))
    iniGrad = 0._r8
end do
```
Now we can set each coordinate with $x_1$%f = 1.23 and so on. Once this is done we can define our function $f(\vec{x})$ and printout the results:
```
type(FG_t)  :: f

f = SQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
f%print()
```
This print for the FG_t type the function value and its gradient (for FGH_t and FGL_t it also prints the Hessian and Lapalcian respectively)

value:

$f$

gradient:

$\left(\frac{f(\vec{x})}{x_1}, \frac{f(\vec{x})}{y_1}, \frac{f(\vec{x})}{z_1}, \frac{f(\vec{x})}{x_2}, \frac{f(\vec{x})}{y_2}, \frac{f(\vec{x})}{z_2}\right)^\text{T}$

If you want to use FGH_t or FGL_t you need to initialize the Hessian to all zeroes or the Laplacian to 0:

FGH_t:
```
iniGrad = 0._r8 !Has size 6
iniHess = 0._r8 !Has size 6x6

do i = 1, SIZE(x)
    iniGrad(i) = 1._r8
    x(i) = FGH_t(0._r8, iniGrad, iniHess) !x has size 6 and represents the positons vector (type(FGH_t) :: x(6))
    iniGrad = 0._r8
end do
```
FGL_t:
```
iniGrad = 0._r8 !Has size 6

do i = 1, SIZE(x)
    iniGrad(i) = 1._r8
    x(i) = FGH_t(0._r8, iniGrad, 0._r8) !x has size 6 and represents the positons vector (type(FGL_t) :: x(6))
    iniGrad = 0._r8
end do
```
There are examples for each type. They will be explained in the next section.

## Examples
3 Examples are included, one for the FG module (value and gradient), one for the FGL module (value, gradient, and Laplacian), and one for the FGH module (value, gradient, and Hessian). The Example code is the optimization of Lennard-Jones clusters and it is commented. If you wish to use the modules for your own code, you will understand how to use it by reading the comments in the examples codes. In order to compile the examples use the following command (compiler (latest version): ifort, ifx, or gfortran; module: FG_m.F90, FGL_m.F90, FGH_m.F90; Example: ExampleFG.f90, ExampleFGL.f90, ExampleFGH.f90; Mathlibrary: -qmkl, -llapack):
```
<compiler> <module> <Example> <Mathlibrary>
```
Examples:
```
ifort FGH_m.F90 ExampleFGH.f90 -qmkl
ifx FGH_m.F90 ExampleFGH.f90 -qmkl
gfortran FGH_m.F90 ExampleFGH.f90 -llapack
```

No other compilers have been tested, so use them at your own risk! 
## Remarks
This program is work in progress! The code has been tested by myself, however there is the possibility of wrong results! If there are any bugs or questions about the code and/or how to use it, feel free to write me.

NEW Release v0.2.0:
1. The FG_t type is now approximately 3 times faster, FGH_t, and FGL_t also are significantly faster!
2. Use of the c preprocessor to shorten the code and to facilitate editing it.
3. Implementation of new operations.
