using AstroLib: trueanom, kepler_solver
using PhysicalConstants
using Unitful, UnitfulAstro
using KeywordCalls
using KeywordDispatch
using Rotations

# Domain specific unit conversions / Constants
const G_nom = 2942.2062175044193 # Rsun^3/Msun/d^2
const G_unit = PhysicalConstants.CODATA2018.G
convert_ρₛ(ρₛ) = ustrip(u"Msun/Rsun^3", ρₛ * u"g/cm^3")

"""
    KeplerianOrbit(; kwargs...)

Keplerian orbit parameterized by the basic observables of a transiting 2-body system.

# Parameters
* `a` - The semi-major axis, nominally in AU
* `aRs`/`aRₛ` - The ratio of the semi-major axis to the star radius.
* `b` - The impact parameter, bounded between 0 ≤ b ≤ 1
* `ecc` - The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1
* `period`/`P` - The orbital period of the planet, nominally in days.
* `rho_star`/`ρₛ` - The spherical star density, nominally in g/cm^3.
* `r_star`/`Rₛ` - The star mass, nominally in solar radii.
* `t0`/`t₀` - The midpoint time of the reference transit, same units as `P`.
* `incl` - The inclination of the orbital plane relative to the axis perpendicular to the
           reference plane, nominally in radians
* `Omega`/`Ω` - The longitude of the ascending node, same units as `incl`.
* `omega`/`ω` - The argument of periapsis, same units as `Ω`.
"""
struct KeplerianOrbit{T,L,D,R,A,I,M} <: AbstractOrbit
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
    Mₛ::M
    aₛ::L
    Mₚ::M
    aₚ::L
end

function normalize_inputs(a, aRₛ, b, ecc, P, Rₛ, t₀, tₚ, t_ref, Mₛ, aₛ, Mₚ, aₚ)
    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # Normalize quantities
    if !(P isa Real)
        a, aₛ, aₚ, Rₛ, = uconvert.(u"Rsun", (a, aₛ, aₚ, Rₛ))
        Mₛ, Mₚ = uconvert.(u"Msun", (Mₛ, Mₚ))
        P, t₀, tₚ, t_ref = uconvert.(u"d", (P, t₀, tₚ, t_ref))
    else
        a, aₛ, aₚ, Rₛ = promote(a, aₛ, aₚ, Rₛ)
        P, t₀, tₚ, t_ref = promote(P, t₀, tₚ, t_ref)
    end

    return a, aRₛ, b, ecc, P, Rₛ, t₀, tₚ, t_ref, Mₛ, aₛ, Mₚ, aₚ
end

function normalize_inputs(
    a::T, aRₛ, b, ecc, P, Rₛ::T,
    t₀, tₚ, t_ref, Mₛ::T, aₛ::T, Mₚ::T, aₚ::T
    ) where {T <: Nothing}

    # Normalize unitless types
    aRₛ, b, ecc = promote(aRₛ, b, ecc)

    # Normalize quantities
    if !(P isa Real)
        P, t₀, tₚ, t_ref = uconvert.(u"d", (P, t₀, tₚ, t_ref))
    else
        P, t₀, tₚ, t_ref = promote(P, t₀, tₚ, t_ref)
    end

    return aRₛ, b, ecc, P, t₀, tₚ, t_ref
end

function _KeplerianOrbit(ρₛ, Rₛ, P, ecc, t₀, incl, Ω, ω)
    # Apply domain specific unit conversions
    ρₛ isa Real && (ρₛ = convert_ρₛ(ρₛ))
    G = P isa Real ? G_nom : G_unit

    # Compute remaining system parameters
    aRₛ = compute_aRₛ(ρₛ, P, G)
    a = compute_a(ρₛ, P, Rₛ, G)
    b = compute_b(ρₛ, P, G, sincos(incl), ecc, ω)
    n = 2.0 * π / P
    M₀ = compute_M₀(ecc, ω)
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀
    Mₛ, aₛ, Mₚ, aₚ = compute_RV_params(ρₛ, Rₛ, a, P, G)

    # Normalize inputs
    a, aRₛ, b, ecc, P, Rₛ, t₀, tₚ, t_ref, Mₛ, aₛ, Mₚ, aₚ = normalize_inputs(
        a, aRₛ, b, ecc, P, Rₛ, t₀, tₚ, t_ref, Mₛ, aₛ, Mₚ, aₚ
    )

    return KeplerianOrbit(a, aRₛ, b, ecc, P, ρₛ, Rₛ, n, t₀, tₚ, t_ref, incl, Ω, ω, M₀,
                          Mₛ, aₛ, Mₚ, aₚ)
end
function _KeplerianOrbit(aRₛ, P, incl, t₀, ecc, Ω, ω)
    # Apply domain specific unit conversions
    G = P isa Real ? G_nom : G_unit

    # Compute remaining system parameters
    ρₛ = compute_ρₛ(aRₛ, P, G)
    b = compute_b(aRₛ, sincos(incl), ecc, ω)
    M₀ = compute_M₀(ecc, ω)
    n = 2.0 * π / P
    tₚ = t₀ - M₀ / n
    t_ref = tₚ - t₀
    a = Rₛ = Mₛ = aₛ = Mₚ = aₚ = nothing

    # Normalize inputs
    aRₛ, b, ecc, P, t₀, tₚ, t_ref = normalize_inputs(
        a, aRₛ, b, ecc, P, Rₛ, t₀, tₚ, t_ref, Mₛ, aₛ, Mₚ, aₚ
    )

    return KeplerianOrbit(a, aRₛ, b, ecc, P, ρₛ, Rₛ, n, t₀, tₚ, t_ref, incl, Ω, ω, M₀,
                          Mₛ, aₛ, Mₚ, aₚ)
end

function KeplerianOrbit(p)
    params = keys(p)
    if :ρₛ ∈ params
        if (:incl ∈ params) & (:b ∉ params)
            return _KeplerianOrbit(p.ρₛ, p.Rₛ, p.P, p.ecc, p.t₀, p.incl, p.Ω, p.ω)
        elseif (:b ∈ params) & (:incl ∉ params)
            G = p.P isa Real ? G_nom : G_unit
            incl = compute_incl(p.ρₛ, p.P, G, p.b, p.ecc, sincos(p.ω))
            return _KeplerianOrbit(p.ρₛ, p.Rₛ, p.P, p.ecc, p.t₀, incl, p.Ω, p.ω)
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end
    elseif :aRₛ ∈ params
        if (:incl ∈ params) & (:b ∉ params)
            return _KeplerianOrbit(p.aRₛ, p.P, p.incl, p.t₀, p.ecc, p.Ω, p.ω)
        elseif (:b ∈ params) & (:incl ∉ params)
            incl = compute_incl(p.aRₛ, p.b, p.ecc, sincos(p.ω))
            return _KeplerianOrbit(p.aRₛ, p.P, incl, p.t₀, p.ecc, p.Ω, p.ω)
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end
    else
        throw(ArgumentError("Please specify either ρₛ or aRₛ"))
    end
end

# KeywordCalls.jl
KeplerianOrbit_KC(nt::NamedTuple{(:ρₛ, :Rₛ, :P, :ecc, :t₀, :incl, :Ω, :ω)}) = _KeplerianOrbit(
    nt.ρₛ, nt.Rₛ, nt.P, nt.ecc, nt.t₀, nt.incl, nt.Ω, nt.ω
)
@kwcall KeplerianOrbit_KC(ρₛ, Rₛ, P, ecc, t₀, incl, Ω, ω)
#@kwalias KeplerianOrbit_KC [
#    rho_s => ρₛ,
#    R_s => Rₛ,
#    period => P,
#    t0 => t₀,
#    Omega => Ω,
#    omega => ω,
#]

# KeywordDispatch.jl
@kwdispatch KeplerianOrbit_KD()
@kwmethod KeplerianOrbit_KD(;ρₛ, Rₛ, P, ecc, t₀, incl, Ω, ω) = _KeplerianOrbit(
    ρₛ, Rₛ, P, ecc, t₀, incl, Ω, ω
)

#############
# Orbit logic
#############
# Star density
compute_ρₛ(aRₛ, P, G_nom) = (3.0 * π / (G_nom * P^2.0)) * aRₛ^3.0
compute_ρₛ(aRₛ, P, G::typeof(G_unit)) = (3.0 * π / (G * P^2.0)) * aRₛ^3.0
compute_ρₛ(a, P, Rₛ, G) = compute_ρₛ(aRₛ(a, Rₛ), P, G)

# Semi-major axis / star radius ratio
compute_aRₛ(a, Rₛ) = a / Rₛ
compute_aRₛ(ρₛ, P, G_nom) = cbrt(G_nom * P^2.0 * ρₛ / (3.0 * π))
compute_aRₛ(ρₛ, P, G::typeof(G_unit)) = cbrt(G * P^2.0 * ρₛ / (3.0 * π))
compute_aRₛ(a, P, Rₛ, G) = aRₛ(compute_ρₛ(a, P, Rₛ, G), P)

# Semi-major axis
compute_a(aRₛ, Rₛ) = aRₛ * Rₛ
compute_a(M_tot, P, G) = cbrt(G * M_tot * P^2 / (4.0 * π^2))
compute_a(ρₛ, P, Rₛ, G) = compute_a(compute_aRₛ(ρₛ, P, G), Rₛ)

# Impact parameter
function compute_b(aRₛ, sincos_incl, ecc, ω)
    sin_ω, cos_ω = sincos(ω)
    incl_factor_inv  = (1.0 - ecc^2.0) / (1.0 + ecc * sin_ω)
    return aRₛ * sincos_incl[2] * incl_factor_inv
end
compute_b(ρₛ, P, G, sincos_incl, ecc, ω) = compute_b(compute_aRₛ(ρₛ, P, G), sincos_incl, ecc, ω)

# Inclination
function compute_incl(aRₛ, b, ecc, sincosω)
    return acos((b/aRₛ) * (1.0 + ecc*sincosω[1])/(1.0 - ecc^2))
end
compute_incl(ρₛ, P, G, b, ecc, sincosω) = compute_incl(compute_aRₛ(ρₛ, P, G), b, ecc, sincosω)

###########
# RV params
###########
compute_M_tot(a, P, G) = 4.0 * π^2.0 * a^3.0 / (G * P^2.0)

function compute_RV_params(ρₛ, Rₛ, a, P, G; Mₚ = zero(typeof(ρₛ * Rₛ^3.0)))
    M_tot = compute_M_tot(a, P, G)
    Mₛ = M_tot - Mₚ
    aₚ = -(Mₛ / M_tot) * a
    aₛ = a + aₚ
    return Mₛ, aₛ, Mₚ, aₚ
end

function compute_M₀(ecc, ω)
    sin_ω, cos_ω = sincos(ω)
    E₀ = 2.0 * atan(√(1.0 - ecc) * cos_ω, √(1.0 + ecc) * (1.0 + sin_ω))
    M₀ = E₀ - ecc * sin(E₀)
    return M₀
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
# separation: aRₛ, aₛ / Rₛ, or aₚ / Rₛ
# TODO: consider moving this to a separate orbital calculations package in the future
function _position(orbit, separation, t)
    sin_ν, cos_ν = compute_true_anomaly(orbit, t)
    if iszero(orbit.ecc)
        r = separation
    else
        r = separation * (1 - orbit.ecc^2) / (1 + orbit.ecc * cos_ν)
    end
    # Transform from orbital plane to equatorial plane
    X = SA[r * cos_ν, r * sin_ν, zero(r)]
    R = RotZXZ(orbit.Ω, -orbit.incl, orbit.ω)
    return R * X
end
_star_position(orb, Rₛ, t) = _position.(orb, orb.aₛ / Rₛ, t)
_planet_position(orb, Rₛ, t) = _position.(orb, orb.aₚ / Rₛ, t)
relative_position(orbit::KeplerianOrbit, t) = _position(orbit, -orbit.aRₛ, t)

# Returns sin(ν), cos(ν)
function compute_true_anomaly(orbit::KeplerianOrbit, t)
    M = (t - orbit.t₀ - orbit.t_ref) * orbit.n
    if iszero(orbit.ecc)
        return sincos(M)
    else
        E = kepler_solver(M, orbit.ecc)
        return sincos(trueanom(E, orbit.ecc))
    end
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

stringify_units(value::Unitful.AbstractQuantity, unit) = value
stringify_units(value, unit) = "$value $unit"
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    print(
        io,
        """
        Keplerian Orbit
          a: $(stringify_units(orbit.a, "R⊙"))
          aRₛ: $(orbit.aRₛ)
          b: $(orbit.b)
          ecc: $(orbit.ecc)
          P: $(stringify_units(orbit.P, "d"))
          ρₛ: $(stringify_units(orbit.ρₛ, "M⊙/R⊙³"))
          Rₛ: $(stringify_units(orbit.Rₛ, "R⊙"))
          t₀: $(stringify_units(orbit.t₀, "d"))
          tₚ: $(stringify_units(orbit.tₚ, "d"))
          t_ref: $(stringify_units(orbit.t_ref, "d"))
          incl: $(stringify_units(orbit.incl, "rad"))
          Ω: $(stringify_units(orbit.Ω, "rad"))
          ω: $(stringify_units(orbit.ω, "rad"))
          Mₛ: $(stringify_units(orbit.Mₛ, "M⊙"))
          aₛ: $(stringify_units(orbit.aₛ, "R⊙"))
          Mₚ: $(stringify_units(orbit.Mₚ, "M⊙"))
          aₚ: $(stringify_units(orbit.aₚ, "R⊙"))
      """
    )
end
