abstract type derivativemethod end
struct backwardsfinite <: derivativemethod end
struct forwardsfinite <: derivativemethod end
struct centralfinite <: derivativemethod end

abstract type nummet end
@kwdef struct ForwardsEuler <:nummet
	name = "Forwards Euler"
	order = 1
	explicit = true
end
@kwdef struct RK4 <: nummet
	name = "RK4"
	order = 4
	explicit = true

end

struct f2Point
	t
	u
end
