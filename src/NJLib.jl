module NJLib

using Base: @kwdef
using CodeTracking

include("types.jl")

include("utils.jl")

include("HighOrderSolve.jl")

include("solver.jl")

include("Numerics.jl")

export Square, Transpose, Det, IdMatrix, Invert
export derivativemethod, backwardsfinite, forwardsfinite, centralfinite
export nummet, ForwardsEuler, RK4
export f2point
export deriv, Newton, fixpoint
export Nstep
export HOProb, initproblem, displayproblem
export problem, solve
end
