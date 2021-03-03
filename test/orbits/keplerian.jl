using PyCall
using Transits.Orbits: KeplerianOrbit, relative_position

pyimport("pip")["main"](["install", "batman-package", "numpy==1.20.1"])

@testset "KeplerianOrbit: sky coords" begin
    # Comparison coords from `batman`
    py"""
    import numpy as np
    from batman import _rsky

    sky_test = {}

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

        sky_test["m_sum"] = m.sum()
        sky_test["r_batman"] = r_batman
        sky_test["m"] = m
        sky_test["t"] = t
        sky_test["t0"] = t0
        sky_test["period"] = period
        sky_test["a"] = a
        sky_test["e"] = e
        sky_test["omega"] = omega
        sky_test["incl"] = incl

        return sky_test
    """
    sky_test = py"sky_coords"()
    allclose = py"np.allclose"

    # Return length(t) × (x; y; z) matrix
    function compute_xyz(orbit, t)
        pos = relative_position.(orbit, t)
        return reduce(hcat, pos)'
    end

    # Create comparison orbits from Transits.jl
    orbits = [
        KeplerianOrbit(
            aRₛ = sky_test["a"][i],
            P = sky_test["period"][i],
            t₀ = sky_test["t0"][i],
            ecc = sky_test["e"][i],
            ω = sky_test["omega"][i],
            incl = sky_test["incl"][i],
        )
        for i in 1:length(sky_test["t0"])
    ]

    # Compute coords
    x = Matrix{Float64}(undef, length(sky_test["t"]), length(sky_test["t0"]))
    y = similar(x)
    z = similar(x)
    for (orbit, x_i, y_i, z_i) in zip(orbits, eachcol(x), eachcol(y), eachcol(z))
        a, b, c = eachcol(compute_xyz(orbit, sky_test["t"]))
        x_i .= a
        y_i .= b
        z_i .= c
    end

    # Compare
    m = sky_test["m"]
    r = @. √(x^2 + y^2)
    r_Transits = r[m]
    r_batman = sky_test["r_batman"][m]

    @test sum(m) > 0
    @test allclose(r_Transits, r_batman, atol=2e-5)
    @test all(z[m] .> 0)
    no_transit = @. (z[!(m)] < 0) | (r[!(m)] > 2)
    @test all(no_transit)
end
