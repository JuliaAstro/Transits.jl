using Conda
using PyCall
using Transits.Orbits: KeplerianOrbit, relative_position
using Unitful, UnitfulAstro

Conda.add(["numpy", "batman-package"]; channel="conda-forge")

# Common Python functionality
py"""
import numpy as np
from batman import _rsky
"""
const allclose = py"np.allclose"

# Common Julia functionality
function compute_r(orbit, t)
    pos = relative_position.(orbit, t)
    r = map(pos) do arr
        √(arr[1]^2 + arr[2]^2)
    end
    return r
end

@testset "KeplerianOrbit: sky coords" begin
    # Comparison coords from `batman`
    py"""
    sky_coords = {}

    def sky_coords():
        t = np.linspace(-100, 100, 1_000)

        t0, period, a, e, omega, incl = (
            x.flatten()
            for x in np.meshgrid(
                np.linspace(-5.0, 5.0, 2),
                np.exp(np.linspace(np.log(5.0), np.log(50.0), 3)),
                np.linspace(50.0, 100.0, 2),
                np.linspace(0.0, 0.9, 5),
                np.linspace(-np.pi, np.pi, 3),
                np.arccos(np.linspace(0, 1, 5)[:-1]),
            )
        )

        r_batman = np.empty((len(t), len(t0)))

        for i in range(len(t0)):
            r_batman[:, i] = _rsky._rsky(
                t, t0[i], period[i], a[i], incl[i], e[i], omega[i], 1, 1
            )

        m = r_batman < 100.0

        return {
            "m_sum" : m.sum(),
            "r_batman" : r_batman,
            "m" : m,
            "t" : t,
            "t0" : t0,
            "period" : period,
            "a" : a,
            "e" : e,
            "omega" : omega,
            "incl" : incl,
        }
    """
    sky_coords = py"sky_coords"()

    # Return length(t) × (x; y; z) matrix
    function compute_xyz(orbit, t)
        pos = relative_position.(orbit, t)
        return reduce(hcat, pos)'
    end

    # Create comparison orbits from Transits.jl
    orbits = [
        KeplerianOrbit(
            aRₛ = sky_coords["a"][i],
            P = sky_coords["period"][i],
            t₀ = sky_coords["t0"][i],
            ecc = sky_coords["e"][i],
            ω = sky_coords["omega"][i],
            incl = sky_coords["incl"][i],
        )
        for i in 1:length(sky_coords["t0"])
    ]

    # Compute coords
    x = Matrix{Float64}(undef, length(sky_coords["t"]), length(sky_coords["t0"]))
    y = similar(x)
    z = similar(x)
    for (orbit, x_i, y_i, z_i) in zip(orbits, eachcol(x), eachcol(y), eachcol(z))
        a, b, c = eachcol(compute_xyz(orbit, sky_coords["t"]))
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

@testset "KeplerianOrbit: small star" begin
    # Model inputs
    m_star = 0.151 # Solar masses
    r_star = 0.189 # Solar radii
    period = 0.4626413 # Days
    t0 = 0.2 # Days
    b = 0.5
    ecc = 0.1
    ω = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        Rₛ = ustrip(u"cm", r_star * u"Rsun"),
        Mₛ = ustrip(u"g", m_star * u"Msun"),
        P =  ustrip(u"s", period * u"d"),
        t₀ = ustrip(u"s", t0 * u"d"),
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
    t_days = test_vals["t"]
    t = ustrip.(u"s", t_days * u"d")
    r_batman = test_vals["r_batman"]
    m = test_vals["m"]
    r = compute_r(orbit, t)
    @test sum(m) > 0
    @test allclose(r_batman[m], r[m], atol=2e-5)
end

@testset "KeplerianOrbit: impact" begin
    # Model inputs
    m_star = 0.151 # Solar masses
    r_star = 0.189 # Solar radii
    period = 0.4626413 # Days
    t0 = 0.2 # Days
    b = 0.5
    ecc = 0.8
    ω = 0.1

    # Sample model from `Transits.jl`
    orbit = KeplerianOrbit(
        Rₛ = ustrip(u"cm", r_star * u"Rsun"),
        Mₛ = ustrip(u"g", m_star * u"Msun"),
        P =  ustrip(u"s", period * u"d"),
        t₀ = ustrip(u"s", t0 * u"d"),
        b = b,
        ecc = ecc,
        ω = ω,
    )

    pos = relative_position.(orbit, orbit.t₀)
    @test allclose((√(pos[1]^2 + pos[2]^2)), orbit.b)
end
