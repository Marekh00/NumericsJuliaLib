@kwdef mutable struct HOProb
	initial::Array{Float64}
	interval
	ODE
	inhom
	order
	step = 0.01
	method = RK4()
	sol_u
	sol_t = Float64[interval[1]]
	probMatrix
end

function initproblem(;initial,
					interval,
					ODE,
					inhomFunc=0,
					step=0.01)
	initial = Float64.(initial)
	order = length(initial)

	solvedUMatrix = Matrix{Float64}(undef,order,1)
	solvedUMatrix[:,1]=initial

	
	problemMatrix = zeros(Float64,order,order)
	for i in 1:order-1
		problemMatrix[i,i+1]=1
	end
	problemMatrix[end,:] = ODE[1:end-1] .* (-1/ODE[end])
	
	prob = HOProb(
		initial = initial,
		interval = interval,
		ODE = ODE,
		inhom = inhomFunc,
		order = order,
		sol_u = solvedUMatrix,
		probMatrix = problemMatrix
	)
	println("Matriserepresentasjon av ODE")
	display(prob.probMatrix)
	return prob
end

function displayproblem(prob::HOProb)
	println("Matrix representation:")
	display(prob.probMatrix)
	if prob.inhom != 0
		println(" = ", @code_string prob.inhom(1))
	end
	println("Initial values: ", prob.initial)
	println("Interval: ", prob.interval)
end
