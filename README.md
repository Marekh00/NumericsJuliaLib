## Numerics-Julia Library

A handcrafted library for solving ODE's

### Features:
- Solver for n-th order linear ODE's with constant coefficients
- Both homogeneous and inhomogeneous

### Method:
The library takes the given ODE and reduces it to a system of
equations which are solved with RK4

### Usage:
Define an ODE inital-value-problem using the `initproblem` function:

KWargs:
- Initial:
Vector of length order-1 with inital values.
The order of initial values goes: [y, y', y'', ... , y^(n-1)]
