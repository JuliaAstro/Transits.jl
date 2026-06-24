using ParallelTestRunner: runtests, find_tests, parse_args
using Transits

const init_code = quote
    using ChainRulesCore
    using ChainRulesTestUtils
    using Orbits
    using QuadGK
    using StableRNGs
    using Test
    using Transits
    using Unitful, UnitfulAstro

    const PLOT = get(ENV, "TEST_PLOTS", "false") == "true"
    const rng = StableRNG(2752)

    PLOT && include("plots.jl")

    @eval Unitful.preferunits(u"Rsun,Msun,d"...)
    ENV["UNITFUL_FANCY_EXPONENTS"] = false

    # Numpy version of `isapprox`
    # https://stackoverflow.com/questions/27098844/allclose-how-to-check-if-two-arrays-are-close-in-julia/27100515#27100515
    function allclose(a, b; rtol=1e-5, atol=1e-8)
        return mapreduce((a, b) -> abs(a - b) ≤ (atol + rtol * abs(b)), &, a, b)
    end
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(Transits, args; testsuite, init_code)
