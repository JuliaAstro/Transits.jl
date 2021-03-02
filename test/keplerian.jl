using Transits.Orbits: relative_position, in_transit

@testset "keplerian orbit: in transits" begin
    t = np.linspace(-20, 20, 1000)
    r_star = 1.5
    orbit = KeplerianOrbit(
        r_star=r_star,
        t0=np.array([0.5, 17.4]),
        period=np.array([10.0, 5.3]),
        ecc=np.array([0.1, 0.8]),
        omega=np.array([0.5, 1.3]),
        Omega=np.array([0.0, 1.0]),
        m_planet=m_planet,
    )

    r_pl = np.array([0.1, 0.03])
    coords = theano.function([], orbit.get_relative_position(t))()
    r2 = coords[0] ** 2 + coords[1] ** 2
    inds = theano.function([], orbit.in_transit(t, r=r_pl))()

    m = np.isin(np.arange(len(t)), inds)
    in_ = r2[inds] <= ((r_star + r_pl) ** 2)[None, :]
    in_ &= coords[2][inds] > 0
    assert np.all(np.any(in_, axis=1))

    out = r2[~m] > ((r_star + r_pl) ** 2)[None, :]
    out |= coords[2][~m] <= 0
    assert np.all(out)
end
