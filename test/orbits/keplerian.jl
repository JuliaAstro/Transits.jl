using BenchmarkTools
using Unitful, UnitfulAstro
using Transits.Orbits: KeplerianOrbit, flip,
                       compute_a, compute_incl,
                       relative_position,
                       _star_position, _planet_position

# https://stackoverflow.com/questions/27098844/allclose-how-to-check-if-two-arrays-are-close-in-julia/27100515#27100515
function allclose(a, b; rtol=1e-5, atol=1e-8)
    return all(abs.(a - b) .<= (atol .+ rtol * abs.(b)))
end

function compute_r(orbit, t)
    pos = relative_position.(orbit, t)
    r = map(pos) do arr
        hypot(arr[1], arr[2])
    end
    return r
end

# Convert vector of vectors -> matrix
as_matrix(pos) = reinterpret(reshape, Float64, pos) |> permutedims

# Tests from:
# https://github.com/exoplanet-dev/exoplanet/blob/main/tests/orbits/keplerian_test.py

@testset "KeplerianOrbit: sky coords" begin
    # Comparison coords from `batman`
    sky_coords = load("./python_code/test_data/KeplerianOrbit_sky_coords.jld2")

    # Create comparison orbits from Transits.jl
    orbits = [
        KeplerianOrbit(
            aR_star = sky_coords["a"][i],
            P = sky_coords["period"][i],
            incl = sky_coords["incl"][i],
            t_0 = sky_coords["t0"][i],
            ecc = sky_coords["e"][i],
            Omega = 0.0,
            omega = sky_coords["omega"][i],
        )
        for i in 1:length(sky_coords["t0"])
    ]

    # Compute coords
    t = sky_coords["t"]
    x = Matrix{Float64}(undef, length(sky_coords["t"]), length(sky_coords["t0"]))
    y = similar(x)
    z = similar(x)
    for (orbit, x_i, y_i, z_i) in zip(orbits, eachcol(x), eachcol(y), eachcol(z))
        pos = relative_position.(orbit, t) |> as_matrix
        a, b, c = eachcol(pos)
        x_i .= a
        y_i .= b
        z_i .= c
    end

    # Compare
    m = sky_coords["m"]
    r = @. √(x^2 + y^2)
    r_Transits = r[m]
    r_batman = sky_coords["r_batman"][m]

    @test sum(m) > 0
    @test allclose(r_Transits, r_batman, atol=2e-5)
    @test all(z[m] .> 0)
    no_transit = @. (z[!(m)] < 0) | (r[!(m)] > 2)
    @test all(no_transit)
end

@testset "KeplerianOrbit: construction performance" begin
    b_ρₛ = @benchmark KeplerianOrbit(
        ρₛ = 2.0,
        R_star = 0.5,
        period = 2.0,
        ecc = 0.0,
        t_0 = 0.0,
        incl = π / 2.0,
        Omega = 0.0,
        omega = 0.0,
    )

    b_ρₛ_units = @benchmark KeplerianOrbit(
        ρₛ = 2.0u"g/cm^3",
        R_star = 0.5u"Rsun",
        period = 2.0u"d",
        ecc = 0.0,
        t_0 = 0.0u"d",
        incl = 90.0u"°",
        Omega = 0.0u"°",
        omega = 0.0u"°",
    )

    b_aR_star = @benchmark KeplerianOrbit(
        aR_star = 7.5,
        P = 2.0,
        incl = π / 2.0,
        t_0 = 0.0,
        ecc = 0.0,
        Omega = 0.0,
        omega = 0.0,
    )

    b_aR_star_units = @benchmark KeplerianOrbit(
        aR_star = 7.5,
        P = 2.0u"d",
        incl = 90.0u"°",
        t_0 = 0.0u"d",
        ecc = 0.0,
        Omega = 0.0u"°",
        omega = 0.0u"°",
    )

    # Units
    @test median(b_ρₛ_units.times) ≤ 100_000 # ns
    @test b_ρₛ_units.allocs ≤ 500
    @test median(b_aR_star_units.times) ≤ 100_000   # ns
    @test b_aR_star_units.allocs ≤ 500

    if v"1.6" ≤ Base.VERSION < v"1.7-"
        @test b_ρₛ.allocs == b_ρₛ.memory == 0
        @test median(b_ρₛ.times) ≤ 500 # ns
        @test b_aR_star.allocs == b_aR_star.memory == 0
        @test median(b_aR_star.times) ≤ 500 # ns
    else
        # TODO: investigate performance regression
        @test median(b_ρₛ.times) ≤ 10_000 # ns
        @test median(b_aR_star.times) ≤ 10_000 # ns
    end
end

@testset "KeplerianOrbit: helper functions" begin
    a, R_star = 2.0, 4.0
    period, G_nom = √π, 1.0
    rho_star = 3.0 * 0.5^3
    b = 0.0
    ecc = 0.0
    sincosomega = (1.0, 0.0)
    @test compute_incl(a, R_star, b, ecc, sincosomega) ≈ π/2.0
end

#=
@testset "KeplerianOrbit: valid inputs" begin
    @test KeplerianOrbit(
        ρₛ = 2.0,
        R_star = 0.5,
        period = 2.0,
        ecc = 0.0,
        t_0 = 0.0,
        incl = π / 2.0,
        Omega = 0.0,
        omega = 0.0,
    ) ===
    KeplerianOrbit(
        ρₛ = 2.0,
        R_star = 0.5,
        period = 2.0,
        ecc = 0.0,
        t_0 = 0.0,
        b = 0.0,
        Omega = 0.0,
        omega = 0.0,
    )
    @test KeplerianOrbit(
        aR_star = 7.5,
        P = 2.0,
        incl = π / 2.0,
        t_0 = 0.0,
        ecc = 0.0,
        Omega = 0.0,
        omega = 0.0,
    ) ===
    KeplerianOrbit(
        aR_star = 7.5,
        P = 2.0,
        b = 0.0,
        t_0 = 0.0,
        ecc = 0.0,
        Omega = 0.0,
        omega = 0.0,
    )
end

@testset "KeplerianOrbit: small star" begin
    # Model inputs
    r_star = 0.189
    m_star = 0.151
    period = 0.4626413
    t0 = 0.2
    b = 0.5
    ecc = 0.1
    omega = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        R_star = r_star,
        Mₛ = m_star,
        P =  period,
        t_0 = t0,
        b = b,
        ecc = ecc,
        omega = omega,
    )

    # Sample model from `batman`
    py"""
    def small_star():
        t = np.linspace(0, $period, 500)
        r_batman = _rsky._rsky(
            t,
            $t0,
            $period,
            $(orbit.aR_star),
            $(orbit.incl),
            $ecc,
            $omega,
            1,
            1
        )

        m = r_batman < 100.0

        return {
            "t": t,
            "r_batman": r_batman,
            "m": m,
        }
    """
    small_star = py"small_star"

    # Compare
    test_vals = small_star()
    t = test_vals["t"]
    r_batman = test_vals["r_batman"]
    m = test_vals["m"]
    r = compute_r(orbit, t)
    @test sum(m) > 0
    @test allclose(r_batman[m], r[m], atol=2e-5)
end

@testset "KeplerianOrbit: impact" begin
    # Model inputs
    r_star = 0.189
    m_star = 0.151
    period = 0.4626413
    t0 = 0.2
    b = 0.5
    ecc = 0.8
    omega = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        R_star = r_star,
        Mₛ = m_star,
        P =  period,
        t_0 = t0,
        b = b,
        ecc = ecc,
        omega = omega,
    )

    pos = relative_position.(orbit, orbit.t_0)
    @test allclose((√(pos[1]^2 + pos[2]^2)), orbit.b)
end
=#

@testset "KeplerianOrbit: flip" begin
    orbit = KeplerianOrbit(
        M_star = 1.3,
        R_star = 1.1,
        t_0 = 0.5,
        period = 100.0,
        ecc = 0.3,
        incl = 0.25*π,
        omega = 0.5,
        Omega = 1.0,
        M_planet = 0.1,
    )

    orbit_flipped = flip(orbit, 0.7)

    t = range(0, 100; length=1_000)

    u_star = as_matrix(_star_position.(orbit, orbit.R_star, t))
    u_planet_flipped = as_matrix(_planet_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end

    u_planet = as_matrix(_planet_position.(orbit, orbit.R_star, t))
    u_star_flipped = as_matrix(_star_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
    end
end

@testset "KeplerianOrbit: flip circular" begin
    t = range(0, 100; length=1_000)

    orbit = KeplerianOrbit(
        M_star = 1.3,
        M_planet = 0.1,
        R_star = 1.0,
        P = 100.0,
        t_0 = 0.5,
        incl = 45.0,
        ecc = 0.0,
        omega = 0.5,
        Omega = 1.0
    )
    orbit_flipped = flip(orbit, 0.7)

    u_star = as_matrix(_star_position.(orbit, orbit.R_star, t))
    u_planet_flipped = as_matrix(_planet_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end

    u_planet = as_matrix(_planet_position.(orbit, orbit.R_star, t))
    u_star_flipped = as_matrix(_star_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
    end
end
