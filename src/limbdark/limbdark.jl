
abstract type AbstractLimbDark end

# struct LimbDarkLightCurve{CT<:AbstractVector,DT}
#     u::CT
#     c::CT
#     c_norm::CT
#     driver::DT
# end

# @doc raw"""
#     LimbDarkLightCurve(u)

# Limb-darkened light curve with limb-darkening coefficients `u`. The limb darkening law can be specThe limb-darkening profile has the quadratic form

# ```math
# \frac{I(\mu)}{I(1)}=1-u_1(1-\mu)-u_2(1-\mu)^2
# ```
# """
# function LimbDarkLightCurve(u::AbstractVector{T}) where T
#     u_ext = pushfirst!(u, -one(T))
#     c = get_cl(u_ext)
#     c_norm = c ./ (π * (c[1] + 2 * c[2] / 3))
#     ld = GreensLimbDark(length(c) - 1)
#     return LimbDarkLightCurve(u, c, c_norm, ld)
# end

(ld::AbstractLimbDark)(args...; kwargs...) = compute(ld, args...; kwargs...)

function compute(ld::AbstractLimbDark, orbit::AbstractOrbit, t, r)
    coords = Orbits.relative_position(orbit, t)
    z = sqrt(coords[1]^2 + coords[2]^2)
    b = z / orbit.r_star
    return compute(ld, b, r)
end

export PolynomialLimbDark, IntegratedLimbDark

include("elliptic.jl")
include("poly.jl")
include("integrated.jl")
