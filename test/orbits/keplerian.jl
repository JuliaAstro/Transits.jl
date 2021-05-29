using BenchmarkTools
using Transits.Orbits: KeplerianOrbit, flip,
                       _star_position, _planet_position, relative_position

function compute_r(orbit, t)
    pos = relative_position.(orbit, t)
    r = map(pos) do arr
        hypot(arr[1], arr[2])
    end
    return r
end

# Tests from:
# https://github.com/exoplanet-dev/exoplanet/blob/main/tests/orbits/keplerian_test.py

@testset "KeplerianOrbit: sky coords" begin
    # Comparison coords from `batman`
    sky_coords = load("./python_code/test_data/KeplerianOrbit_sky_coords.jld2")

    # Convert vector of vectors -> matrix
    as_matrix(pos) = reinterpret(reshape, Float64, pos) |> permutedims

    # Create comparison orbits from Transits.jl
    orbits = [
        KeplerianOrbit(
            aRₛ = sky_coords["a"][i],
            P = sky_coords["period"][i],
            incl = sky_coords["incl"][i],
            t₀ = sky_coords["t0"][i],
            ecc = sky_coords["e"][i],
            Ω = 0.0,
            ω = sky_coords["omega"][i],
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

@testset "KeplerianOrbit: construction" begin
    b_ρₛ = @benchmark KeplerianOrbit(
        ρₛ = 2.0,
        Rₛ = 0.5,
        period = 2.0,
        ecc = 0.0,
        t₀ = 0.0,
        incl = π / 2.0,
        Ω = 0.0,
        ω = 0.0,
    )
    b_aRₛ = @benchmark KeplerianOrbit(
        aRₛ = 7.5,
        P = 2.0,
        incl = π / 2.0,
        t₀ = 0.0,
        ecc = 0.0,
        Ω = 0.0,
        ω = 0.0,
    )

    @test b_ρₛ.allocs == b_ρₛ.memory == 0
    @test median(b_ρₛ.times) ≤ 500 # ns
    @test b_aRₛ.allocs == b_aRₛ.memory == 0
    @test median(b_aRₛ.times) ≤ 500 # ns
end

#=
@testset "KeplerianOrbit: small star" begin
    # Model inputs
    r_star = 0.189
    m_star = 0.151
    period = 0.4626413
    t0 = 0.2
    b = 0.5
    ecc = 0.1
    ω = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        Rₛ = r_star,
        Mₛ = m_star,
        P =  period,
        t₀ = t0,
        b = b,
        ecc = ecc,
        ω = ω,
    )

    # Sample model from `batman`
    py"""
    def small_star():
        t = np.linspace(0, $period, 500)
        r_batman = _rsky._rsky(
            t,
            $t0,
            $period,
            $(orbit.aRₛ),
            $(orbit.incl),
            $ecc,
            $ω,
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
    ω = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        Rₛ = r_star,
        Mₛ = m_star,
        P =  period,
        t₀ = t0,
        b = b,
        ecc = ecc,
        ω = ω,
    )

    pos = relative_position.(orbit, orbit.t₀)
    @test allclose((√(pos[1]^2 + pos[2]^2)), orbit.b)
end

@testset "KeplerianOrbit: flip" begin
    stack(arr_arr) = hcat((reshape(map(p -> p[i], arr_arr), :) for i in 1:3)...)

    t = range(0, 100; length=1_000)

    orbit = KeplerianOrbit(
        Mₛ = 1.3,
        Mₚ = 0.1,
        Rₛ = 1.0,
        P = 100.0,
        t₀ = 0.5,
        incl = 45.0,
        ecc = 0.3,
        ω = 0.5,
        Ω = 1.0
    )
    orbit_flipped = flip(orbit, 0.7)


    u_star = stack(_star_position.(orbit, orbit.Rₛ, t))
    u_planet_flipped = stack(_planet_position.(orbit_flipped, orbit.Rₛ, t))
    for i in 1:3
        @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end

    u_planet = stack(_planet_position.(orbit, orbit.Rₛ, t))
    u_star_flipped = stack(_star_position.(orbit_flipped, orbit.Rₛ, t))
    for i in 1:3
        @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
    end
end

@testset "KeplerianOrbit: flip circular" begin
    stack(arr_arr) = hcat((reshape(map(p -> p[i], arr_arr), :) for i in 1:3)...)

    t = range(0, 100; length=1_000)

    orbit = KeplerianOrbit(
        Mₛ = 1.3,
        Mₚ = 0.1,
        Rₛ = 1.0,
        P = 100.0,
        t₀ = 0.5,
        incl = 45.0,
        ecc = 0.0,
        ω = 0.5,
        Ω = 1.0
    )
    orbit_flipped = flip(orbit, 0.7)

    u_star = stack(_star_position.(orbit, orbit.Rₛ, t))
    u_planet_flipped = stack(_planet_position.(orbit_flipped, orbit.Rₛ, t))
    for i in 1:3
        @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end

    u_planet = stack(_planet_position.(orbit, orbit.Rₛ, t))
    u_star_flipped = stack(_star_position.(orbit_flipped, orbit.Rₛ, t))
    for i in 1:3
        @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
    end
end
=#
