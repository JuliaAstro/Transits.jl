module Orbits

using StaticArrays

abstract type AbstractOrbit end


relative_position(::AbstractOrbit, t)


include("simple.jl")

end # module
