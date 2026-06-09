@kwdef struct problem
	initial::Array{Float64}
	interval
	func
	step = 0.01
	method = RK4()
	sol_t = Float64[interval[1]]
	sol_u = Float64[initial[1]]
end


function solve(prob::problem)
	nsteps = round(Int, (prob.interval[2]-prob.interval[1])/prob.step)
	for i in 1:nsteps
		a,b = Nstep(prob.method; func=prob.func, tn=prob.sol_t[end], un=prob.sol_u[end], step=prob.step)
		push!(prob.sol_t, a)
		push!(prob.sol_u, b)
	end
end

function solve(prob::HOProb)
	nsteps = round(Int, (prob.interval[2]-prob.interval[1])/prob.step)
	C = zeros(prob.order)
	C[end]=1/prob.ODE[end]
	un = prob.initial[:]

	new_sol_t = Vector{Float64}(undef, nsteps+1)
	new_sol_t[1] = prob.sol_t[1]

	new_sol_u = similar(prob.sol_u, prob.order, nsteps+1)
	new_sol_u[:, 1] = prob.sol_u[:, 1]

	if prob.inhom == 0
		for i in 1:nsteps
			tnp, un = Nstep(
			prob.method, Val{:homogeneous}(); mat = prob.probMatrix,
			tn = new_sol_t[i],
			un = un,
			step = prob.step)
			new_sol_t[i+1] = tnp
			new_sol_u[:, i+1] = un
		end
	else
		for i in 1:nsteps
			tnp, un = Nstep(
			prob.method, Val{:system}(); mat = prob.probMatrix,
			inhomfunc = prob.inhom,
			adderterm = C,
			tn = new_sol_t[i],
			un = un,
			step = prob.step)
			new_sol_t[i+1] = tnp
			new_sol_u[:, i+1] = un
		end
	end
	prob.sol_t = new_sol_t
	prob.sol_u = new_sol_u
end
