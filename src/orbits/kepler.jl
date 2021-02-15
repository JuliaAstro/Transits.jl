
using AstroLib: kepler_solver, trueanom

struct KeplerianOrbit <: AbstractOrbit
    a
    ecc
    period
    ρ_star
    r_star
    m_star
    m_planet
    m_total
    n
    a_star
    a_planet
    K0
    M0
    t_periastron
    t0
    tref
    incl
    sin_incl
    cos_incl
    Ω
    sinΩ
    cosΩ
    ω
    sinω
    cosω

end

const G_grav = 2942.2062175044193

"""
    KeplerOrbit(; period, duration, t0=0, b=0, r_star=1)
Keplerian orbit parameterized by the basic observables of a transiting 2-body system.
# Parameters
* `period` - The orbital period of the planets, nominally in days
* `m_planet` -  The planet mass, nominally in Jupiter masses
* `r_star` - The star mass, nominally in solar radii
* `m_star` - The star mass, nominally in solar masses
* `ρ_star` - The spherical star density, nominally in g/cc
* `ecc` - The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1
* `t_periastron` - The time of periastron, same units as `period`
* `t0` - The midpoint time of the reference transit, same units as `period`
* `Ω` - The longitude of the ascending node, nominally in radians
* `ω` - The argument of periapsis, same units as `Ω`
"""
function KeplerianOrbit(;
    period,
    m_planet=0,
    r_star=1,
    m_star=1,
    ρ_star = 3 * m_star / (4 * π * r_star^3),
    ecc=nothing,
    t_periastron=0,
    t0 = 0 + t_periastron,
    incl = 0,
    Ω = 3*π / 2,
    ω = π / 2,
)
    a = cbrt(G_grav * (m_star + m_planet) * period^2 / (4 * π^2))
    m_total = m_star + m_planet
    n = 2 * π / period
    a_star = a * m_planet / m_total
    a_planet = -a * m_star / m_total
    K0 = n * a / m_total

    M0 = 0.5 * π
    if ecc === nothing # Circular orbit (e == 0)
        ecc = 0.0
    end

    # If b given, can compute i from ω
    # ecc_factor = (1. + ecc*sin(ω))/(1. - ecc^2)
    # inc_inv_factor = (b/a)*ecc_factor
    # incl = acos(inc_inv_factor)

    # Euler angles
    sin_incl, cos_incl = sincos(incl)
    sinΩ, cosΩ = sincos(Ω)
    sinω, cosω = sincos(ω)

    tref = t_periastron - t0

    return KeplerianOrbit(
        a,
        ecc,
        period,
        ρ_star,
        r_star,
        m_star,
        m_planet,
        m_total,
        n,
        a_star,
        a_planet,
        K0,
        M0,
        t_periastron,
        t0,
        tref,
        incl,
        sin_incl,
        cos_incl,
        Ω,
        sinΩ,
        cosΩ,
        ω,
        sinω,
        cosω,
    )
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
function relative_position(orbit::KeplerianOrbit, t)
    sinν, cosν = get_true_anomaly(orbit, t)
    if orbit.ecc === nothing
        r = -orbit.a
    else
        r = -orbit.a * (1 - orbit.ecc^2) / (1 + orbit.ecc * cosν)
    end
    return rotate_vector(orbit, r * cosν, r * sinν)
end

# Returns sin(ν), cos(ν)
function get_true_anomaly(orbit::KeplerianOrbit, t)
    M = orbit.n * ((t - orbit.t0) - orbit.tref)
    E = kepler_solver(M, orbit.ecc)
    return sincos(trueanom(E, orbit.ecc))
end
#(M, ::Nothing) = sincos(M)

# Transform from orbital plane to equatorial plane
function rotate_vector(orbit::KeplerianOrbit, x, y)
    # rotate about z0 axis by ω
    if orbit.ecc === nothing
        x1, y1 = x, y
    else
        #=
        # https://en.wikipedia.org/wiki/Orbital_elements
        mat_ω = [
            orbit.cosω orbit.sinω 0
           -orbit.sinω orbit.cosω 0
                 0          0     1
        ]

        mat_incl = [
            1       0              0
            0  orbit.cos_incl orbit.sin_incl
            0 -orbit.sin_incl orbit.cos_incl
        ]

        mat_Ω = [
            orbit.cosΩ orbit.sinΩ 0
           -orbit.sinΩ orbit.cosΩ 0
                 0          0     1
        ]
        =#
        x1 = orbit.cosω * x - orbit.sinω * y
        y1 = orbit.sinω * x + orbit.cosω * y
    end

    # rotate about x1 axis by -incl
    x2 = x1
    y2 = orbit.cos_incl * y1
    Z = -orbit.sin_incl * y1

    # rotate about z2 axis by Ω
     if orbit.Ω === nothing
         return SA[x2, y2, Z]
     end
     X = orbit.cosΩ * x2 - orbit.sinΩ * y2
     Y = orbit.sinΩ * x2 + orbit.cosΩ * y2
     return SA[X, Y, Z]
end
