module Orbits

using StaticArrays

export SimpleOrbit, KeplerianOrbit

abstract type AbstractOrbit end

Base.broadcastable(o::AbstractOrbit) = Ref(o)

"""
    relative_position(::AbstractOrbit, t)
"""
relative_position(::AbstractOrbit, t)


include("simple.jl")
include("kepler.jl")

end # module
