module Transits

using ChainRulesCore
using Orbits
using Orbits: AbstractOrbit
using StaticArrays
using Unitful

export AbstractLimbDark,
       PolynomialLimbDark,
       QuadLimbDark,
       IntegratedLimbDark,
       SecondaryLimbDark,
       compute,
       # distributions
       Kipping13

"""
    AbstractLimbDark

A limb dark law need only need to implement `compute(::Law, b, r)` to extend the limb darkening interface.

# See also
[`compute`](@ref)
"""
abstract type AbstractLimbDark end

Base.broadcastable(law::AbstractLimbDark) = Ref(law)

"""
    (::AbstractLimbDark)(b, r)

An alias for calling [`compute`](@ref)

# Examples
```jldoctest
julia> ld = PolynomialLimbDark([0.4, 0.26]);

julia> ld(0, 0.01)
0.9998785437247428
```
"""
(ld::AbstractLimbDark)(args...; kwargs...) = compute(ld, args...; kwargs...)

"""
    compute(::AbstractLimbDark, b, r; kwargs...)

Compute the relative flux for the given impact parameter `b` and radius ratio `r`. The impact parameter is unitless. The radius ratio is given in terms of the host; e.g., if the companion is half the size of the host, r=0.5.
"""
function compute end

"""
    compute(::AbstractLimbDark, orbit::AbstractOrbit, t, r)

Compute the relative flux by calculating the impact parameter at time `t` from the given orbit. The time needs to be compatible with the period of the orbit, nominally in days.

# Examples
```jldoctest orb
julia> ld = PolynomialLimbDark([0.4, 0.26]);

julia> orbit = SimpleOrbit(period=3, duration=1);

julia> ld(orbit, 0, 0.1) # primary egress
0.9878664434953113

julia> ld(orbit, 0.1, 0.1) # 0.1 d
0.9879670695533511
```
this works effortlessly with libraries like [Unitful.jl](https://github.com/painterqubits/Unitful.jl)

```jldoctest orb
julia> using Unitful

julia> orbit = SimpleOrbit(period=3u"d", duration=3u"hr");

julia> ld(orbit, 0u"d", 0.1)
0.9878664434953113
```
"""
function compute(ld::AbstractLimbDark, orbit::AbstractOrbit, t, r)
    coords = Orbits.relative_position(orbit, t)
    # strip units, if necessary (e.g., KeplerianOrbit)
    x, y, z = ustrip.(coords)
    # make sure los is in front of star
    if z > 0
        μ = sqrt(x^2 + y^2)
        return compute(ld, μ, r)
    else
        return one(Base.promote_typeof(x, y, z))
    end
end

include("polynomial/elliptic.jl")
include("polynomial/series.jl")
include("polynomial/poly.jl")
include("polynomial/poly-grads.jl")
include("polynomial/quad.jl")
include("polynomial/quad-grads.jl")
include("integrated.jl")
include("secondary.jl")
include("distributions.jl")

end
