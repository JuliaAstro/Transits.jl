using Conda
using PyCall
using QuadGK
using StableRNGs
using Transits
using Test
using Unitful, UnitfulAstro

# Use Julia's own version of Python
ENV["PYTHON"] = ""
Conda.add(["numpy", "batman-package"]; channel="conda-forge")

const PLOT = get(ENV, "TEST_PLOTS", "false") == "true"
PLOT && include("plots.jl")

rng = StableRNG(2752)

@testset "Transits" begin
    include("Mn_integral.jl")
    include("distributions.jl")
    include("elliptic.jl")
    include("orbits/keplerian.jl")
    include("poly.jl")
    include("show.jl")
end
