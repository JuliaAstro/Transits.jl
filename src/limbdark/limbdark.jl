
abstract type AbstractLimbDark end

Base.broadcastable(law::AbstractLimbDark) = Ref(law)

(ld::AbstractLimbDark)(args...; kwargs...) = compute(ld, args...; kwargs...)

function compute(ld::AbstractLimbDark, orbit::AbstractOrbit, t, r)
    coords = Orbits.relative_position(orbit, t)
    z = sqrt(coords[1]^2 + coords[2]^2)
    b = z / orbit.r_star
    return compute(ld, b, r)
end

include("elliptic.jl")
include("poly.jl")
include("integrated.jl")
