using PyCall
using QuadGK
using StableRNGs
using Transits
using Test
#using Unitful, UnitfulAstro

py"""
import numpy as np
from batman import _rsky
"""

# Numpy version of `isapprox`
allclose(a, b; rtol=1e-5, atol=1e-8) = all(@. abs(a - b) ≤ (atol + rtol*abs(b)))

const PLOT = get(ENV, "TEST_PLOTS", "false") == "true"
PLOT && include("plots.jl")

rng = StableRNG(2752)

@testset "Transits" begin
    include("Mn_integral.jl")
    include("distributions.jl")
    include("elliptic.jl")
    include("orbits/keplerian.jl")
    include("orbits/solvers.jl")
    include("poly.jl")
    include("show.jl")
end
