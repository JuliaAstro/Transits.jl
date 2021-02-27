using AstroLib: kepler_solver, trueanom
using KeywordDispatch
using PhysicalConstants
using Unitful
using UnitfulAstro

const G = PhysicalConstants.CODATA2018.G
const G_cgs = ustrip(u"cm^3/g/s^2", G)
const M₀ = π / 2.0

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
struct KeplerianOrbit{T,L,D,R,A,I}
    a::L
    aRₛ::R
    b::R
    ecc::R
    P::T
    ρₛ::D
    Rₛ::L
    n::I
    t₀::T
    incl::A
    Ω::R
    ω::R
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
)

@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, ecc, P, t₀, incl)
    # Apply relevant conversions to CGS
    #Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    #P isa Real && (P = P * 86_400.0)
    #incl isa Real && (incl = incl * π / 180.0)
    Ω = M₀
    ω = 0.0
    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
    a = compute_a(ρₛ, P, Rₛ)
    b = compute_b(ρₛ, P, sincos(incl))
    n = 2.0 * π / P

    # Normalize quantities
    a, Rₛ = promote(a, Rₛ)

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

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
        incl,
        Ω,
        ω,
    )
end

@kwmethod function KeplerianOrbit(;aRₛ, P, b, t₀, ecc)
    Ω = M₀
    ω = 0.0
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))

    return KeplerianOrbit(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        nothing,
        nothing,
        2.0 * π / P,
        t₀,
        incl,
        Ω,
        ω,
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
@kwmethod compute_aRₛ(;ρₛ, P::T) where {T <: Unitful.Time} = cbrt(G * P^2.0 * ρₛ / (3.0 * π))
@kwmethod compute_aRₛ(;a, P, Rₛ) = aRₛ(compute_ρₛ(a, P, Rₛ), P)
@kwmethod compute_aRₛ(;a, Rₛ) = a / Rₛ

# Semi-major axis
compute_a(ρₛ, P, Rₛ) = compute_a(compute_aRₛ(ρₛ=ρₛ, P=P), Rₛ)
compute_a(aRₛ, Rₛ) = aRₛ * Rₛ

# Impact parameter
compute_b(ρₛ, P, sincosi) = compute_b(compute_aRₛ(ρₛ=ρₛ, P=P), sincosi)
compute_b(aRₛ, sincosi) = aRₛ * sincosi[2]

# Inclination
function compute_incl(aRₛ, b, ecc, sincosω)
    return acos((b/aRₛ) * (1.0 + ecc*sincosω[1])/(1.0 - ecc^2))
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
function relative_position(orbit::KeplerianOrbit, t)
    sinν, cosν = compute_true_anomaly(orbit, t)
    if iszero(orbit.ecc)
        r = orbit.aRₛ
    else
        r = orbit.aRₛ * (1 - orbit.ecc^2) / (1 + orbit.ecc * cosν)
    end
    return rotate_vector(orbit, r * cosν, r * sinν)
end

function compute_M(orbit::KeplerianOrbit, t)
    if iszero(orbit.ecc)
        M = M₀ + orbit.n * (t - orbit.t₀) / 2.0
    else
        sinω, cosω = sincos(orbit.ω)
        E₀ = 2.0 * atan(
            sqrt(1.0 - orbit.ecc) * cosω,
            sqrt(1.0 + orbit.ecc) * (1.0 + sinω),
        )
        M₀_ecc = E₀ - orbit.ecc * sin(E₀)
        M = M₀_ecc + orbit.n * (t - orbit.t₀)
    end
    return M
end

# Returns sin(ν), cos(ν)
function compute_true_anomaly(orbit::KeplerianOrbit, t)
    #M = compute_M(orbit, t)
    M = M₀ + orbit.n * (t - orbit.t₀) / 2.0
    E = kepler_solver(M, orbit.ecc)
    return sincos(trueanom(E, orbit.ecc))
end

# Transform from orbital plane to equatorial plane
function rotate_vector(orbit::KeplerianOrbit, x, y)
    sini, cosi = sincos(orbit.incl)
    sinΩ, cosΩ = sincos(orbit.Ω)
    sinω, cosω = sincos(orbit.ω)

    # Rotate about z0 axis by ω
    if iszero(orbit.ecc)
        x1, y1 = x, y
    else
        x1 = cosω * x - sinω * y
        y1 = sinω * x + cosω * y
    end

    # Rotate about x1 axis by -incl
    x2 = x1
    y2 = cosi * y1
    Z = -sini * y1

    # Rotate about z2 axis by Ω
    X = cosΩ * x2 - sinΩ * y2
    Y = sinΩ * x2 + cosΩ * y2

    return SA[X, Y, Z]
end

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
    a = orbit.a isa Quantity ? uconvert(u"AU", orbit.a) : orbit.a
    print(
        io,
        """
        KeplerianOrbit
         a: $a
         aRₛ: $(upreferred(orbit.aRₛ))
         b: $(upreferred(orbit.b))
         ecc: $(orbit.ecc)
         P: $(orbit.P)
         ρₛ: $(orbit.ρₛ)
         Rₛ: $(orbit.Rₛ)
         t₀: $(orbit.t₀)
         incl: $(orbit.incl)
         Ω: $(orbit.Ω)
         ω: $(orbit.ω)
        """
    )
end
