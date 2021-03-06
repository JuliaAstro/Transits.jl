using Transits.Orbits: trueanom, kepler_solver

# Common Python functionality
py"""
import numpy as np
"""
const allclose = py"np.allclose"

# Tests from:
# https://github.com/dfm/kepler.py/blob/main/tests/test_kepler.py
@testset "kepler_solver: M, E edge case" begin
    # Compute sin_ν, cos_ν without using arctan function directly
    function compute_sincos_ν_no_atan(E, ecc)
        sin_E, cos_E = sincos(E)
        denom = 1.0 + cos_E
        x = √((1 + ecc) / (1 - ecc)) * sin_E / denom # where ν = 2 arctan(x)
        x² = x * x
        # Apply trig identites:
        #     sincos(arctan x) = (x, 1) ./ √(1 + x²)
        #     sincos(ν) = 2 sin(arctan x)cos(arctan x), 1 - 2 sin²(arctan x)
        denom = 1.0 / (1.0 + x²)
        sin_ν = 2.0 * x * denom
        cos_ν = (1.0 - x²) * denom
        return sin_ν , cos_ν
    end

    function compute_E_solver(E, ecc)
        M =  E - ecc * sin(E)
        return kepler_solver(M, ecc) # <- E
    end

    function edge(E, ecc)
        # Restrict E to the range [-π, π]
        E = rem2pi(E, RoundNearest)
        # Computed from solver
        E_solver = compute_E_solver(E, ecc)
        println(E_solver)
        sin_ν_solver, cos_ν_solver = compute_sincos_ν_no_atan(E_solver, ecc)
        # Computed directly
        sin_ν, cos_ν = sincos(trueanom(E, ecc))
        return [E_solver, E, sin_ν_solver, sin_ν, cos_ν_solver, cos_ν]
    end

    # Test different inputs for E and ecc
    Es = [0.0, 2 * π, -226.2, -170.4]
    eccs = fill(1.0 - 1e-6, length(Es))
    eccs[end] = 0.9939879759519037
    (
        E_solver,     E_user,
        sin_ν_solver, sin_ν_user,
        cos_ν_solver, cos_ν_user,
    ) = eachrow(hcat(edge.(Es, eccs)...));

    @test all(isfinite.(sin_ν_solver))
    @test all(isfinite.(cos_ν_solver))
    @test allclose(E_solver, E_user)
    @test allclose(sin_ν_solver, sin_ν_user)
    @test allclose(cos_ν_solver, cos_ν_user)
end
