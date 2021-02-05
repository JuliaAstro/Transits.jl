module Orbits

using StaticArrays

export SimpleOrbit

abstract type AbstractOrbit end

Base.broadcastable(o::AbstractOrbit) = Ref(o)

"""
    relative_position(::AbstractOrbit, t)
"""
relative_position(::AbstractOrbit, t)


include("simple.jl")

end # module
