using Transits.Orbits: trueanom, kepler_solver

# Compute sin_ν, cos_ν without using arctan function directly
function compute_sincos_ν_no_atan(E, ecc; tol=1e-10)
    sin_E, cos_E = sincos(E)
    denom = 1.0 + cos_E

    # Adjust denominator if necessary to avoid dividing by zero
    m = denom > tol
    m_inv = !m
    denom += 1.0 * m_inv

    # Compute sin ν, cos ν
    x = √((1 + ecc) / (1 - ecc)) * sin_E / denom # where ν = 2 arctan(x)
    x² = x * x
    # Apply trig identites:
    #     sincos(arctan x) = (x, 1) ./ √(1 + x²)
    #     sincos(ν) = 2 sin(arctan x)cos(arctan x), 1 - 2 sin²(arctan x)
    denom = 1.0 / (1.0 + x²)
    sin_ν = 2.0 * x * denom
    cos_ν = (1.0 - x²) * denom
    return sin_ν * m , cos_ν * m - 1.0 * m_inv
end

function compute_E_solver(E, ecc)
    M =  E - ecc * sin(E)
    return kepler_solver(M, ecc) # <- E
end

function _compute_vals(E, ecc)
    # Restrict E to the range [-π, π]
    E = rem2pi(Float64(E), RoundNearest)
    # Computed from solver
    E_solver = compute_E_solver(E, ecc)
    sin_ν_solver, cos_ν_solver = compute_sincos_ν_no_atan(E_solver, ecc)
    # Computed directly
    sin_ν, cos_ν = sincos(trueanom(E, ecc))
    return [E_solver, E, sin_ν_solver, sin_ν, cos_ν_solver, cos_ν]
end
# Returns a generator that can be unpacked into:
#     (E_solver,    E_user,
#     sin_ν_solver, sin_ν_user,
#     cos_ν_solver, cos_ν_user)
compute_summary(Es, eccs) = eachrow(mapreduce(_compute_vals, hcat, Es, eccs))

function test_vals(summary)
    (E_solver,     E_user,
     sin_ν_solver, sin_ν_user,
     cos_ν_solver, cos_ν_user) = summary

    @test all(isfinite.(sin_ν_solver))
    @test all(isfinite.(cos_ν_solver))
    @test allclose(E_solver, E_user)
    @test allclose(sin_ν_solver, sin_ν_user)
    @test allclose(cos_ν_solver, cos_ν_user)
end

# Tests from:
# https://github.com/dfm/kepler.py/blob/main/tests/test_kepler.py

@testset "kepler_solver: M, E edge case" begin
    Es = [0.0, 2 * π, -226.2, -170.4]
    eccs = fill(1.0 - 1e-6, length(Es))
    eccs[end] = 0.9939879759519037
    test_vals(compute_summary(Es, eccs))
end

@testset "kepler_solver: Let them eat π" begin
    eccs = range(0.0, 1.0; length=100)[begin:end-1]
    Es = fill(π, length(eccs))
    test_vals(compute_summary(Es, eccs))
end

@testset "kepler_solver: solver" begin
    ecc_range = range(0, 1; length=500)[begin:end-1]
    E_range = range(-300, 300; length=1_001)
    E_ecc_pairs = Iterators.product(E_range, ecc_range)
    Es = reshape(map(x -> x[1], E_ecc_pairs), :, 1)
    eccs = reshape(map(x -> x[2], E_ecc_pairs), :, 1)
    test_vals(compute_summary(Es, eccs))
end
