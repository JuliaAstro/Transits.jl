
# using LinearAlgebra: dot

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


# function (lc::LimbDarkLightCurve)(orbit::AbstractOrbit, t, r)

#     coords = Orbits.relative_position(orbit, t)
#     b = sqrt(coords[1]^2 + coords[2]^2)
#     if coords[3] > 0
#         b_ = abs(b) / orbit.r_star
#         r_ = abs(r) / orbit.r_star
#         if b_ < 1 + r_
#             sT = lc.driver(b_, r_)
#             return dot(sT, lc.c_norm) - 1
#         end
#     end
#     return zero(b)
# end

export PolynomialLimbDark

include("elliptic.jl")
include("poly.jl")
# include("integrate.jl")
