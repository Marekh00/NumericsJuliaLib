## Numerics-Julia Library

A handcrafted library for solving ODE's

### Features:
- Solver for n-th order linear ODE's with constant coefficients
- Both homogeneous and inhomogeneous
- fuction displayproblem(problem) to display problem-data (will be improved on later)

### Method:
The library takes the given ODE and reduces it to a system of
equations which are solved with RK4

### Usage:
Define an ODE inital-value-problem using the `initproblem` function:
Gives an HOProb-object. Usefull variables contained in object are:
- probMatrix:
Matrixrepresentation of ODE-problem

- sol_t:
x-axis or time-values used in solving

- sol_u:
vector containing solved y, y', y'' ..., y^(n-1)

KWargs:
- Initial:
Vector of length order-1 with inital values.
The order of initial values goes: [y, y', y'', ... , y^(n-1)]

- Interval:
Tuple with startpoint and endpoint values. Assumed that inital values are from startpoint.

- ODE:
Vector of coefficients of ODE.
Coefficients are in order [y, y', y''. ... , y^n]

Optional:
inhomFunc:
Known function of the ODE in terms of t


### Example:
y''' + y'' + y' + y = x^2
```julia
using NJLib
initialV = [1,2,1]
interval = (0,5)
ODE = [1,1,1,1]
f(t) = t^2
prob = initproblem(initial = initalV, interval = interval, ODE = ODE, inhomFunc = f)
solve(prob)
```
Solution is now contained in prob.sol_u.
(specifically, y-values are in prob.sol_u[1,:])

### NOTE:
function solve and problem-struct is for the moment deprecated. It's use is mainly for 1.order ODE, and will be brought up to standard later.
