
using AstroLib: kepler_solver

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
end

const G_grav = 2942.2062175044193

function KeplerianOrbit(; period, m_planet=0, r_star=1, m_star=1, ρ_star = 3 * m_star / (4 * π * r_star^3), ecc=nothing, t_periastron=0, t0 = 0 + t_periastron)
    a = cbrt(G_grav * (m_star + m_planet) * period^2 / (4 * π^2))
    m_total = m_star + m_planet
    n = 2 * π / period
    a_star = a * m_planet / m_total
    a_planet = -a * m_star / m_total
    K0 = n * a / m_total
    if ecc === nothing
        M0 = 0.5 * π
        incl_factor = 1
    else
        error("not implemented")
    end
    incl = 0.5 * π
    sin_incl, cos_incl = sincos(incl)

    tref = t_periastron - t0
    return KeplerianOrbit(a,
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
                          cos_incl)
end

function relative_position(orbit::KeplerianOrbit, t)
    sinf, cosf = get_true_anomaly(orbit, t)
    if orbit.ecc === nothing
        r = -orbit.a
    else
        r = -orbit.a * (1 - orbit.ecc^2) / (1 + orbit.ecc * cosf)
    end
    return rotate_vector(orbit, r * cosf, r * sinf)
end

function get_true_anomaly(orbit::KeplerianOrbit, t)
    M = orbit.n * ((t - orbit.t0) - orbit.tref)
    return kepler(M, orbit.ecc)
end

function kepler(M, ecc)
    E = kepler_solver(M, ecc)
    return sincos(E)
end
kepler(M, ::Nothing) = sincos(M)

function rotate_vector(orbit::KeplerianOrbit, x, y)
    # rotate about z0 axis by ω
    if orbit.ecc === nothing
        x1, y1 = x, y
    else
        # x1 = orbit.cosω * x - orbit.sinω * y
        # y1 = orbit.sinω * x + orbit.cosω * y
    end

    # rotate about x1 axis by -incl
    x2 = x1
    y2 = orbit.cos_incl * y1
    Z = -orbit.sin_incl * y1

    return SA[x2, y2, Z]

    # rotate about z2 axis by Ω
    # if orbit.Ω === nothing
    #     return SA[x2, y2, Z]
    # end
    # X = orbit.cosΩ * x2 - orbit.sinΩ * y2
    # Y = orbit.sinΩ * x2 + orbit.cosΩ * y2
    # return SA[X, Y, Z]
end