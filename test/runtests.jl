using ChainRulesCore
using ChainRulesTestUtils
using Orbits
using QuadGK
using StableRNGs
using Test
using Transits
using Unitful, UnitfulAstro

Unitful.preferunits(u"Rsun,Msun,d"...)
ENV["UNITFUL_FANCY_EXPONENTS"] = false


# Numpy version of `isapprox`
# https://stackoverflow.com/questions/27098844/allclose-how-to-check-if-two-arrays-are-close-in-julia/27100515#27100515
allclose(a, b; rtol=1e-5, atol=1e-8) = mapreduce((a,b) -> abs(a - b) ≤ (atol + rtol*abs(b)), &, a, b)

const PLOT = get(ENV, "TEST_PLOTS", "false") == "true"
PLOT && include("plots.jl")

rng = StableRNG(2752)

@testset "Transits" begin
    include("Mn_integral.jl")
    include("distributions.jl")
    include("elliptic.jl")
    include("grads.jl")
    include("poly.jl")
    include("show.jl")
end
