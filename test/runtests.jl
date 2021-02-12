using StableRNGs
using Transits
using Test

rng = StableRNG(2752)

include("elliptic.jl")
include("poly_transit.jl")
