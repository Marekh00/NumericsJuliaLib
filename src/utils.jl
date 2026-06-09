#
# Her er div. viktige funksjoner for numeriske metoder
#


function deriv(method::backwardsfinite, func, x; step = 0.01)
	return (func(x) - func(x-step))/step
end

function deriv(method::forwardsfinite, func, x; step = 0.01)
	return (func(x+step)-func(x))/step
end

function deriv(method::centralfinite, func, x; step = 0.01)
	return (func(x+step)-func(x-step))/(2*step)
end

function Newton(func, x0; step=0.01, iter::Int = 20)
	xn = 0
	for i in 1:iter
		xn = x0 - func(x0)/deriv(backwardsfinite(),func,x0,step=step)
		x0 = xn
	end
	return x0
end

function fixpoint(func, x0; iter::Int = 20)
	xn = 0
	for i in 1:iter
		xn = func(x0)
		x0 = xn
	end
	return x0
end

function Nstep(m::ForwardsEuler;func, tn, un, step = 0.01)
	un = un + step*func(tn,un)
	return tn+step,un
end

function Nstep(m::RK4;func, tn, un, step = 0.01)
	K1 = func(tn,			un)
	K2 = func(tn+0.5*step,	un+0.5*step*K1)
	K3 = func(tn+0.5*step,	un+0.5*step*K2)
	K4 = func(tn+step,		un+step*K3)
	un += (step/6)*(K1+2*K2+2*K3+K4)
	return tn+step,un
end

function Nstep(m::RK4, ::Val{:system}; mat, inhomfunc, adderterm, tn, un, step=0.01)
	
	K1 = mat*(un)+adderterm*inhomfunc(tn)
	K2 = mat*(un+0.5*step*K1)+adderterm*inhomfunc(tn+0.5*step)
	K3 = mat*(un+0.5*step*K2) + adderterm*inhomfunc(tn+0.5*step)
	K4 = mat*(un+step*K3) + adderterm*inhomfunc(tn+step)
	un += (step/6)*(K1+2*K2+2*K3+K4)
	return tn+step,un
end

function Nstep(m::RK4, ::Val{:homogeneous}; mat, tn, un, step=0.01)
	K1 = mat*un
	K2 = mat*(un+0.5*step*K1)
	K3 = mat*(un+0.5*step*K2)
	K4 = mat*(un+step*K3)
	un += (step/6)*(K1+2*K2+2*K3+K4)
	return tn+step,un
end

