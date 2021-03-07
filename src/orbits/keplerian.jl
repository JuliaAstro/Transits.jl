using AstroLib: trueanom, kepler_solver
using KeywordDispatch
using PhysicalConstants
using Unitful
using UnitfulAstro

const G = PhysicalConstants.CODATA2018.G
const G_cgs = ustrip(u"cm^3/g/s^2", G)

"""
    KeplerianOrbit(; kwargs...)
Keplerian orbit parameterized by the basic observables of a transiting 2-body system.
# Parameters
* `a` - The semi-major axis, nominally in AU
* `aRₛ` - The ratio of the semi-major axis to the star radius. Aliased to `aRs`
* `b` - The impact parameter, bounded between 0 ≤ b ≤ 1
* `ecc` - The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1
* `P` - The orbital period of the planet, nominally in days
* `ρₛ` - The spherical star density, nominally in g/cc. Aliased to `rho_s`
* `Rₛ` - The star mass, nominally in solar radii. Aliased to `R_s`
* `t₀` - The midpoint time of the reference transit, same units as `P`. Aliased to `t0`
* `incl` - The inclination of the orbital plane relative to the axis perpendicular to the
           reference plane, nominally in degrees
* `Ω` - The longitude of the ascending node, nominally in radians. Aliased to `Omega`
* `ω` - The argument of periapsis, same units as `Ω`. Aliased to `omega`
"""
struct KeplerianOrbit{T,L,D,R,A,I} <: AbstractOrbit
    a::L
    aRₛ::R
    b::R
    ecc::R
    P::T
    ρₛ::D
    Rₛ::L
    n::I
    t₀::T
    tₚ::T
    t_ref::T
    incl::A
    Ω::R
    ω::R
    M₀::R
    Mₛ
    aₛ
    Mₚ
    aₚ
end

# Enable keyword dispatch and argument name aliasing
@kwdispatch KeplerianOrbit(;
    Omega => Ω,
    omega => ω,
    aRs => aRₛ,
    rho_s => ρₛ,
    aRs => aRₛ,
    Rs => Rₛ,
    t0 => t₀,
    M0 => M₀,
    tp => tₚ,
    Mp => Mₚ,
    Rp => Rₚ,
)

@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, ecc, P, t₀, incl)
    # Apply relevant conversions to CGS
    #Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    #P isa Real && (P = P * 86_400.0)
    #incl isa Real && (incl = incl * π / 180.0)
    Ω = 0.0
    ω = 0.0
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(ρₛ=ρₛ, P=P, Rₛ=Rₛ)
    b = compute_b(ρₛ, P, sincos(incl), ecc, ω)
    n = 2.0 * π / P
    M₀ = compute_M₀(ecc, ω)
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    # Normalize quantities
    a, Rₛ = promote(a, Rₛ)

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, ecc, P, tₚ, incl)
    # Apply relevant conversions to CGS
    #Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    #P isa Real && (P = P * 86_400.0)
    #incl isa Real && (incl = incl * π / 180.0)
    Ω = 0.0
    ω = 0.0
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(ρₛ=ρₛ, P=P, Rₛ=Rₛ)
    b = compute_b(ρₛ, P, sincos(incl), ecc, ω)
    n = 2.0 * π / P
    M₀ = compute_M₀(ecc, ω)
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀

    # Normalize quantities
    a, Rₛ = promote(a, Rₛ)

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, P, t₀, b, ecc, ω)
    # Apply relevant conversions to CGS
    #Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    #P isa Real && (P = P * 86_400.0)
    #incl isa Real && (incl = incl * π / 180.0)
    Ω = 0.0
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(ρₛ=ρₛ, P=P, Rₛ=Rₛ)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    n = 2.0 * π / P
    M₀ = compute_M₀(ecc, ω)
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀
    Mₛ = compute_Mₛ(ρₛ, Mₛ)

    # Normalize quantities
    a, Rₛ = promote(a, Rₛ)

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, P, tₚ, b, ecc, ω)
    # Apply relevant conversions to CGS
    #Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    #P isa Real && (P = P * 86_400.0)
    #incl isa Real && (incl = incl * π / 180.0)
    Ω = 0.0
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(ρₛ=ρₛ, P=P, Rₛ=Rₛ)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    n = 2.0 * π / P
    M₀ = compute_M₀(ecc, ω)
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀
    Mₛ = compute_Mₛ(ρₛ, Mₛ)

    # Normalize quantities
    a, Rₛ = promote(a, Rₛ)

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;aRₛ, P, b, t₀, ecc)
    Ω = 0.0
    ω = 0.0
    ρₛ = compute_ρₛ(aRₛ, P)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    return KeplerianOrbit(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        nothing,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        nothing,
        nothing,
        nothing,
        nothing
    )
end

@kwmethod function KeplerianOrbit(;aRₛ, P, b, tₚ, ecc)
    Ω = 0.0
    ω = 0.0
    ρₛ = compute_ρₛ(aRₛ, P)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀

    return KeplerianOrbit(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        nothing,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        nothing,
        nothing,
        nothing,
        nothing
    )
end

@kwmethod function KeplerianOrbit(;aRₛ, P, t₀, ecc, ω, incl)
    Ω = 0.0
    ρₛ = compute_ρₛ(aRₛ, P)
    b = compute_b(aRₛ, sincos(incl), ecc, ω)
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    return KeplerianOrbit(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        nothing,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        nothing,
        nothing,
        nothing,
        nothing
    )
end

@kwmethod function KeplerianOrbit(;aRₛ, P, tₚ, ecc, ω, incl)
    Ω = 0.0
    ρₛ = compute_ρₛ(aRₛ, P)
    b = compute_b(aRₛ, sincos(incl), ecc, ω)
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀

    return KeplerianOrbit(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        nothing,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        nothing,
        nothing,
        nothing,
        nothing
    )
end

@kwmethod function KeplerianOrbit(;Mₛ, Rₛ, P, t₀, b, ecc, ω)
    Ω = 0.0
    ρₛ = 3.0 * Mₛ / ( 4.0 * π * Rₛ^3.0 )
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(aRₛ=aRₛ, Rₛ=Rₛ)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;Mₛ, Rₛ, P, tₚ, b, ecc, ω)
    Ω = 0.0
    ρₛ = 3.0 * Mₛ / ( 4.0 * π * Rₛ^3.0 )
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(aRₛ=aRₛ, Rₛ=Rₛ)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    t₀ =  tₚ + M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;Mₛ, Mₚ, Rₛ, P, t₀, incl, ecc, ω, Ω)
    ρₛ = 3.0 * Mₛ / ( 4.0 * π * Rₛ^3.0 )
    M_tot = Mₛ + Mₚ
    a = compute_a(M_tot=M_tot, P=P)
    aRₛ = compute_aRₛ(a=a, Rₛ=Rₛ)
    b = compute_b(aRₛ, sincos(incl), ecc, ω)
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P; Mₚ=Mₚ)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;Mₛ, Mₚ, Rₛ, P, tₚ, incl, ecc, ω, Ω)
    ρₛ = 3.0 * Mₛ / ( 4.0 * π * Rₛ^3.0 )
    M_tot = Mₛ + Mₚ
    a = compute_a(M_tot=M_tot, P=P)
    aRₛ = compute_aRₛ(a=a, Rₛ=Rₛ)
    b = compute_b(aRₛ, sincos(incl), ecc, ω)
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P; Mₚ=Mₚ)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;a, Rₛ, P, t₀, b, ecc, ω)
    Ω = 0.0
    aRₛ = compute_aRₛ(a=a, Rₛ=Rₛ)
    ρₛ = compute_ρₛ(aRₛ, P)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

@kwmethod function KeplerianOrbit(;a, Rₛ, P, tₚ, b, ecc, ω)
    Ω = 0.0
    aRₛ = compute_aRₛ(a=a, Rₛ=Rₛ)
    ρₛ = compute_ρₛ(aRₛ, P)
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    t₀ = tₚ + M₀ / n
    t_ref = tₚ - t₀

    # RV info
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P)

    return KeplerianOrbit(
        a,
        aRₛ,
        b,
        ecc,
        P,
        ρₛ,
        Rₛ,
        n,
        t₀,
        tₚ,
        t_ref,
        incl,
        Ω,
        ω,
        M₀,
        Mₛ,
        aₛ,
        Mₚ,
        aₚ
    )
end

#############
# Orbit logic
#############
# Star density
compute_ρₛ(aRₛ, P) = (3.0 * π / (G_cgs * P^2.0)) * aRₛ^3.0
compute_ρₛ(aRₛ, P::T) where {T <: Unitful.Time} = (3.0 * π / (G * P^2.0)) * aRₛ^3.0
compute_ρₛ(a, P, Rₛ) = compute_ρₛ(aRₛ(a, Rₛ), P)

# Semi-major axis / star radius ratio
@kwdispatch compute_aRₛ()
@kwmethod compute_aRₛ(;ρₛ, P) = cbrt(G_cgs * P^2.0 * ρₛ / (3.0 * π))
@kwmethod compute_aRₛ(;ρₛ, P::T) where {T <: Unitful.Time} = cbrt(G_cgs * P^2.0 * ρₛ / (3.0 * π))
@kwmethod compute_aRₛ(;a, P, Rₛ) = aRₛ(compute_ρₛ(a, P, Rₛ), P)
@kwmethod compute_aRₛ(;a, Rₛ) = a / Rₛ

# Semi-major axis
@kwdispatch compute_a()
@kwmethod compute_a(;ρₛ, P, Rₛ) = compute_a(aRₛ=compute_aRₛ(ρₛ=ρₛ, P=P), Rₛ=Rₛ)
@kwmethod compute_a(;aRₛ, Rₛ) = aRₛ * Rₛ
@kwmethod compute_a(;M_tot, P::T) where {T <: Unitful.Time} = cbrt(G * M_tot * P^2 / (4.0 * π^2))
@kwmethod compute_a(;M_tot, P) = cbrt(G_cgs * M_tot * P^2 / (4.0 * π^2))

# Impact parameter
compute_b(ρₛ, P, sincos_incl, ecc, ω) = compute_b(compute_aRₛ(ρₛ=ρₛ, P=P), sincos_incl, ecc, ω)
function compute_b(aRₛ, sincos_incl, ecc, ω)
    sin_ω, cos_ω = sincos(ω)
    incl_factor_inv  = (1.0 - ecc^2.0) / (1.0 + ecc * sin_ω)
    return aRₛ * sincos_incl[2] * incl_factor_inv
end

# Inclination
function compute_incl(aRₛ, b, ecc, sincosω)
    return acos((b/aRₛ) * (1.0 + ecc*sincosω[1])/(1.0 - ecc^2))
end

# RV params
compute_M_tot(a, P) = 4.0 * π^2.0 * a^3.0 / (G_cgs * P^2.0)
compute_M_tot(a, P::T) where {T <: Unitful.Time} = 4.0 * π^2 * a^3 / (G * P^2.0)
function compute_RV_params(ρₛ, Rₛ, a, P; Mₚ = zero(typeof(ρₛ * Rₛ^3.0)))
    M_tot = compute_M_tot(a, P)
    Mₛ = M_tot - Mₚ
    aₚ = -(Mₛ / M_tot) * a
    aₛ = a + aₚ
    return Mₛ, aₛ, Mₚ, aₚ
end

function compute_M₀(ecc, ω)
    if false #iszero(ecc)
        M₀ = π / 2.0
    else
        sinω, cosω = sincos(ω)
        E₀ = 2.0 * atan(
            sqrt(1.0 - ecc) * cosω,
            sqrt(1.0 + ecc) * (1.0 + sinω),
        )
        M₀ = E₀ - ecc * sin(E₀)
    end
    return M₀
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
# a_rel: aRₛ, aₛ / Rₛ, or aₚ / Rₛ
function _position(orbit, a_rel, t)
    sin_ν, cos_ν = compute_true_anomaly(orbit, t)
    if false #iszero(orbit.ecc)
        r = a_rel
    else
        r = a_rel * (1 - orbit.ecc^2) / (1 + orbit.ecc * cos_ν)
    end
    return rotate_vector(orbit, r * cos_ν, r * sin_ν)
end
star_position(orb, Rₛ, t) = _position.(orb, orb.aₛ / Rₛ, t)
planet_position(orb, Rₛ, t) = _position.(orb, orb.aₚ / Rₛ, t)
relative_position(orbit::KeplerianOrbit, t) = _position(orbit, -orbit.aRₛ, t)

# Returns sin(ν), cos(ν)
function compute_true_anomaly(orbit::KeplerianOrbit, t)
    M = orbit.n * (t - orbit.t₀ - orbit.t_ref)
    E = kepler_solver(M, orbit.ecc)
    return sincos(trueanom(E, orbit.ecc))
end

# Transform from orbital plane to equatorial plane
function rotate_vector(orbit::KeplerianOrbit, x, y)
    sin_incl, cos_incl = sincos(orbit.incl)
    sin_Ω, cos_Ω = sincos(orbit.Ω)
    sin_ω, cos_ω = sincos(orbit.ω)

    # Rotate about z0 axis by ω
    if false #iszero(orbit.ecc)
        x1, y1 = x, y
    else
        x1 = cos_ω * x - sin_ω * y
        y1 = sin_ω * x + cos_ω * y
    end

    # Rotate about x1 axis by -incl
    x2 = x1
    y2 = cos_incl * y1
    Z = -sin_incl * y1

    # Rotate about z2 axis by Ω
    X = cos_Ω * x2 - sin_Ω * y2
    Y = sin_Ω * x2 + cos_Ω * y2

    return SA[X, Y, Z]
end

flip(orbit::KeplerianOrbit, Rₚ) = KeplerianOrbit(
    Mₛ = orbit.Mₚ,
    Mₚ = orbit.Mₛ,
    Rₛ = Rₚ,
    P = orbit.P,
    tₚ = orbit.tₚ,
    incl = orbit.incl,
    ecc = orbit.ecc,
    ω = orbit.ω - π,
    Ω = orbit.Ω
)

function Base.show(io::IO, orbit::KeplerianOrbit)
    a = orbit.a isa Nothing ? nothing : uconvert(u"AU", orbit.a)
    aRₛ = orbit.aRₛ
    b = orbit.b
    ecc = orbit.ecc
    P = orbit.P
    ρₛ = orbit.ρₛ
    Rₛ = orbit.Rₛ
    t₀ = orbit.t₀
    incl = orbit.incl
    Ω = orbit.Ω
    ω = orbit.ω
    print(
        io,
        """KeplerianOrbit(
            a=$(upreferred(orbit.a)), aRₛ=$(orbit.aRₛ),
            b=$(orbit.b), ecc=$(orbit.ecc), P=$(orbit.P),
            ρₛ=$(orbit.ρₛ), Rₛ=$(orbit.Rₛ),
            t₀=$(orbit.t₀), incl=$(orbit.incl),
            Ω=$(orbit.Ω), ω = $(orbit.ω)
        )"""
    )
end

function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    a = orbit.a isa Quantity ? orbit.a : "$(ustrip(u"AU", orbit.a * u"cm")) AU"
    P = orbit.P isa Quantity ? orbit.P : "$(ustrip(u"d", orbit.P * u"s")) d"
    ρₛ = orbit.ρₛ isa Quantity ? orbit.ρₛ : "$(orbit.ρₛ) g/cm^3"
    Rₛ = orbit.Rₛ isa Quantity ? orbit.Rₛ : "$(ustrip(u"Rsun", orbit.Rₛ * u"cm")) Rsun"
    t₀ = orbit.t₀ isa Quantity ? orbit.t₀ : "$(ustrip(u"d", orbit.t₀ * u"s")) d"
    tₚ = orbit.tₚ isa Quantity ? orbit.tₚ : "$(ustrip(u"d", orbit.tₚ * u"s")) d"
    t_ref = orbit.t_ref isa Quantity ? orbit.t_ref : "$(ustrip(u"d", orbit.t_ref * u"s")) d"
    incl = orbit.incl isa Quantity ? orbit.incl : "$(ustrip(u"°", orbit.incl))°"
    Ω = orbit.Ω isa Quantity ? orbit.Ω : "$(orbit.Ω) rad"
    ω = orbit.ω isa Quantity ? orbit.ω : "$(orbit.ω) rad"
    Mₛ = orbit.Mₛ isa Quantity ? orbit.Mₛ : "$(ustrip(u"Msun", orbit.Mₛ * u"g")) Msun"
    aₛ = orbit.aₛ isa Quantity ? orbit.aₛ : "$(ustrip(u"AU", orbit.aₛ * u"cm")) AU"
    Mₚ = orbit.Mₚ isa Quantity ? orbit.Mₚ : "$(ustrip(u"Mjup", orbit.Mₚ * u"g")) Mjup"
    aₚ = orbit.aₚ isa Quantity ? orbit.aₚ : "$(ustrip(u"AU", orbit.aₚ * u"cm")) AU"
    print(
        io,
        """
        KeplerianOrbit
         a: $a
         aRₛ: $(upreferred(orbit.aRₛ))
         b: $(upreferred(orbit.b))
         ecc: $(orbit.ecc)
         P: $P
         ρₛ: $ρₛ
         Rₛ: $Rₛ
         t₀: $t₀
         tₚ: $tₚ
         t_ref: $t_ref
         incl: $incl
         Ω: $Ω
         ω: $ω
         Mₛ: $Mₛ
         aₛ: $aₛ
         Mₚ: $Mₚ
         aₚ: $aₚ
        """
    )
end
