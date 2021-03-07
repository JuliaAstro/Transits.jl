using Conda
using PyCall
using QuadGK
using StableRNGs
using Transits
using Test
using Unitful, UnitfulAstro

Conda.add(["numpy", "batman-package"]; channel="conda-forge")
py"""
import numpy as np
from batman import _rsky
"""
const allclose = py"np.allclose"

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
