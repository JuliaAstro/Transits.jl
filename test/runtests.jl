using QuadGK
using StableRNGs
using Transits
using Test

rng = StableRNG(2752)

@testset "Transits" begin
    include("elliptic.jl")
    include("Mn_integral.jl")
    include("poly.jl")
end
