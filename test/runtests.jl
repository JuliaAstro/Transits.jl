using QuadGK
using StableRNGs
using Transits
using Test

rng = StableRNG(2752)

include("elliptic.jl")
include("Mn_integral.jl")
include("poly_transit.jl")
