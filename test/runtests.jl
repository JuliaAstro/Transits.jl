using ChainRulesCore
using ChainRulesTestUtils
using QuadGK
using StableRNGs
using Transits
using Test

const PLOT = get(ENV, "TEST_PLOTS", "false") == "true"
PLOT && include("plots.jl")

rng = StableRNG(2752)

@testset "Transits" begin
    # include("elliptic.jl")
    # include("Mn_integral.jl")
    # include("poly.jl")
    # include("show.jl")
    # include("distributions.jl")
    include("grads.jl")
end
