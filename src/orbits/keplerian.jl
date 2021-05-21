using AstroLib: trueanom, kepler_solver
using PhysicalConstants
using Unitful, UnitfulAstro
using KeywordCalls
using KeywordDispatch
using Rotations

# Domain specific unit conversions / Constants
const G_nom = 2942.2062175044193 # Rsun^3/Msun/d^2
const G_unit = PhysicalConstants.CODATA2018.G
convert_rho_s(rho_s) = ustrip(u"Msun/Rsun^3", rho_s * u"g/cm^3")

"""
    KeplerianOrbit(; kwargs...)

Keplerian orbit parameterized by the basic observables of a transiting 2-body system.

# Parameters
* `a` - The semi-major axis, nominally in AU
* `aR_s`/`aRₛ` - The ratio of the semi-major axis to the star radius.
* `b` - The impact parameter, bounded between 0 ≤ b ≤ 1
* `ecc` - The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1
* `period`/`P` - The orbital period of the planet, nominally in days.
* `rho_s`/`ρₛ` - The spherical star density, nominally in g/cm^3.
* `r_star`/`R_s` - The star mass, nominally in solar radii.
* `t_0`/`t₀` - The midpoint time of the reference transit, same units as `P`.
* `incl` - The inclination of the orbital plane relative to the axis perpendicular to the
           reference plane, nominally in radians
* `Omega`/`Ω` - The longitude of the ascending node, same units as `incl`.
* `omega`/`ω` - The argument of periapsis, same units as `Omega`.
"""
struct KeplerianOrbit{T,L,D,R,A,I,M} <: AbstractOrbit
    a::L
    aR_s::R
    b::R
    ecc::R
    period::T
    rho_s::D
    R_s::L
    n::I
    t_0::T
    t_p::T
    t_ref::T
    incl::A
    Omega::R
    omega::R
    M₀::R
    M_s::M
    a_s::L
    M_p::M
    a_p::L
end

function normalize_inputs(a, aR_s, b, ecc, period, R_s, t_0, t_p, t_ref, M_s, a_s, M_p, a_p)
    # Normalize unitless types
    aR_s, b, ecc = promote(aR_s, b, ecc)

    # Normalize quantities
    if !(period isa Real)
        a, a_s, a_p, R_s, = uconvert.(u"Rsun", (a, a_s, a_p, R_s))
        M_s, M_p = uconvert.(u"Msun", (M_s, M_p))
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
    else
        a, a_s, a_p, R_s = promote(a, a_s, a_p, R_s)
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
    end

    return a, aR_s, b, ecc, period, R_s, t_0, t_p, t_ref, M_s, a_s, M_p, a_p
end

function normalize_inputs(
    a::T, aR_s, b, ecc, period, R_s::T,
    t_0, t_p, t_ref, M_s::T, a_s::T, M_p::T, a_p::T
    ) where {T <: Nothing}

    # Normalize unitless types
    aR_s, b, ecc = promote(aR_s, b, ecc)

    # Normalize quantities
    if !(period isa Real)
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
    else
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
    end

    return aR_s, b, ecc, period, t_0, t_p, t_ref
end

function KeplerianOrbit(rho_s, R_s, period, ecc, t_0, incl, Omega, omega)
    # Apply domain specific unit conversions
    rho_s isa Real && (rho_s = convert_rho_s(rho_s))
    G = period isa Real ? G_nom : G_unit

    # Compute remaining system parameters
    aR_s = compute_aR_s(rho_s, period, G)
    a = compute_a(rho_s, period, R_s, G)
    b = compute_b(rho_s, period, G, sincos(incl), ecc, omega)
    n = 2.0 * π / period
    M₀ = compute_M₀(ecc, omega)
    t_p = t_0 - M₀ / n
    t_ref = t_p - t_0
    M_s, a_s, M_p, a_p = compute_RV_params(rho_s, R_s, a, period, G)

    # Normalize inputs
    a, aR_s, b, ecc, period, R_s, t_0, t_p, t_ref, M_s, a_s, M_p, a_p = normalize_inputs(
        a, aR_s, b, ecc, period, R_s, t_0, t_p, t_ref, M_s, a_s, M_p, a_p
    )

    return KeplerianOrbit(a, aR_s, b, ecc, period, rho_s, R_s, n, t_0, t_p, t_ref, incl, Omega, omega, M₀,
                          M_s, a_s, M_p, a_p)
end

function KeplerianOrbit(aR_s, period, incl, t_0, ecc, Omega, omega)
    # Apply domain specific unit conversions
    G = period isa Real ? G_nom : G_unit

    # Compute remaining system parameters
    rho_s = compute_rho_s(aR_s, period, G)
    b = compute_b(aR_s, sincos(incl), ecc, omega)
    M₀ = compute_M₀(ecc, omega)
    n = 2.0 * π / period
    t_p = t_0 - M₀ / n
    t_ref = t_p - t_0
    a = R_s = M_s = a_s = M_p = a_p = nothing

    # Normalize inputs
    aR_s, b, ecc, period, t_0, t_p, t_ref = normalize_inputs(
        a, aR_s, b, ecc, period, R_s, t_0, t_p, t_ref, M_s, a_s, M_p, a_p
    )

    return KeplerianOrbit(a, aR_s, b, ecc, period, rho_s, R_s, n, t_0, t_p, t_ref, incl, Omega, omega, M₀,
                          M_s, a_s, M_p, a_p)
end

function KeplerianOrbit(nt::NamedTuple{(:rho_s, :R_s, :period, :ecc, :t_0, :incl, :Omega, :omega)})
    params = keys(nt)
    if :rho_s ∈ params
        if (:incl ∈ params) & (:b ∉ params)
            return KeplerianOrbit(nt.rho_s, nt.R_s, nt.period, nt.ecc, nt.t_0, nt.incl, nt.Omega, nt.omega)
        elseif (:b ∈ params) & (:incl ∉ params)
            G = nt.period isa Real ? G_nom : G_unit
            incl = compute_incl(nt.rho_s, nt.period, G, nt.b, nt.ecc, sincos(nt.omega))
            return KeplerianOrbit(nt.rho_s, nt.R_s, nt.period, nt.ecc, nt.t_0, incl, nt.Omega, nt.omega)
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end
    elseif :aR_s ∈ params
        if (:incl ∈ params) & (:b ∉ params)
            return KeplerianOrbit(nt.aR_s, nt.period, nt.incl, nt.t_0, nt.ecc, nt.Omega, nt.omega)
        elseif (:b ∈ params) & (:incl ∉ params)
            incl = compute_incl(nt.aR_s, nt.b, nt.ecc, sincos(nt.omega))
            return KeplerianOrbit(nt.aR_s, nt.period, incl, nt.t_0, nt.ecc, nt.Omega, nt.omega)
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end
    else
        throw(ArgumentError("periodlease specify either rho_s or aR_s"))
    end
end

@kwcall KeplerianOrbit(rho_s, R_s, period, ecc, t_0, incl, Omega, omega)
@kwalias KeplerianOrbit [
    ρₛ => rho_s,
    Rₛ => R_s,
    P => period,
    t₀ => t_0,
    Ω => Omega,
    ω => omega,
]

#############
# Orbit logic
#############
# Star density
compute_rho_s(aR_s, period, G_nom) = (3.0 * π / (G_nom * period^2.0)) * aR_s^3.0
compute_rho_s(aR_s, period, G::typeof(G_unit)) = (3.0 * π / (G * period^2.0)) * aR_s^3.0
compute_rho_s(a, period, R_s, G) = compute_rho_s(aR_s(a, R_s), period, G)

# Semi-major axis / star radius ratio
compute_aR_s(a, R_s) = a / R_s
compute_aR_s(rho_s, period, G_nom) = cbrt(G_nom * period^2.0 * rho_s / (3.0 * π))
compute_aR_s(rho_s, period, G::typeof(G_unit)) = cbrt(G * period^2.0 * rho_s / (3.0 * π))
compute_aR_s(a, period, R_s, G) = aR_s(compute_rho_s(a, period, R_s, G), period)

# Semi-major axis
compute_a(aR_s, R_s) = aR_s * R_s
compute_a(M_tot, period, G) = cbrt(G * M_tot * period^2 / (4.0 * π^2))
compute_a(rho_s, period, R_s, G) = compute_a(compute_aR_s(rho_s, period, G), R_s)

# Impact parameter
function compute_b(aR_s, sincos_incl, ecc, omega)
    sin_omega, cos_omega = sincos(omega)
    incl_factor_inv  = (1.0 - ecc^2.0) / (1.0 + ecc * sin_omega)
    return aR_s * sincos_incl[2] * incl_factor_inv
end
compute_b(rho_s, period, G, sincos_incl, ecc, omega) = compute_b(compute_aR_s(rho_s, period, G), sincos_incl, ecc, omega)

# Inclination
function compute_incl(aR_s, b, ecc, sincosomega)
    return acos((b/aR_s) * (1.0 + ecc*sincosomega[1])/(1.0 - ecc^2))
end
compute_incl(rho_s, period, G, b, ecc, sincosomega) = compute_incl(compute_aR_s(rho_s, period, G), b, ecc, sincosomega)

###########
# RV params
###########
compute_M_tot(a, period, G) = 4.0 * π^2.0 * a^3.0 / (G * period^2.0)

function compute_RV_params(rho_s, R_s, a, period, G; M_p = zero(typeof(rho_s * R_s^3.0)))
    M_tot = compute_M_tot(a, period, G)
    M_s = M_tot - M_p
    a_p = -(M_s / M_tot) * a
    a_s = a + a_p
    return M_s, a_s, M_p, a_p
end

function compute_M₀(ecc, omega)
    sin_omega, cos_omega = sincos(omega)
    E₀ = 2.0 * atan(√(1.0 - ecc) * cos_omega, √(1.0 + ecc) * (1.0 + sin_omega))
    M₀ = E₀ - ecc * sin(E₀)
    return M₀
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
# separation: aR_s, a_s / R_s, or a_p / R_s
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
    R = RotZXZ(orbit.Omega, -orbit.incl, orbit.omega)
    return R * X
end
_star_position(orb, R_s, t) = _position.(orb, orb.a_s / R_s, t)
_planet_position(orb, R_s, t) = _position.(orb, orb.a_p / R_s, t)
relative_position(orbit::KeplerianOrbit, t) = _position(orbit, -orbit.aR_s, t)

# Returns sin(ν), cos(ν)
function compute_true_anomaly(orbit::KeplerianOrbit, t)
    M = (t - orbit.t_0 - orbit.t_ref) * orbit.n
    if iszero(orbit.ecc)
        return sincos(M)
    else
        E = kepler_solver(M, orbit.ecc)
        return sincos(trueanom(E, orbit.ecc))
    end
end

flip(orbit::KeplerianOrbit, Rₚ) = KeplerianOrbit(
    M_s = orbit.M_p,
    M_p = orbit.M_s,
    R_s = Rₚ,
    period = orbit.period,
    t_p = orbit.t_p,
    incl = orbit.incl,
    ecc = orbit.ecc,
    omega = orbit.omega - π,
    Omega = orbit.Omega
)

stringify_units(value::Unitful.AbstractQuantity, unit) = value
stringify_units(value, unit) = "$value $unit"
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    print(
        io,
        """
        Keplerian Orbit
          a: $(stringify_units(orbit.a, "R⊙"))
          aRₛ: $(orbit.aR_s)
          b: $(orbit.b)
          ecc: $(orbit.ecc)
          P: $(stringify_units(orbit.period, "d"))
          ρₛ: $(stringify_units(orbit.rho_s, "M⊙/R⊙³"))
          Rₛ: $(stringify_units(orbit.R_s, "R⊙"))
          t₀: $(stringify_units(orbit.t_0, "d"))
          tₚ: $(stringify_units(orbit.t_p, "d"))
          t_ref: $(stringify_units(orbit.t_ref, "d"))
          incl: $(stringify_units(orbit.incl, "rad"))
          Ω: $(stringify_units(orbit.Omega, "rad"))
          ω: $(stringify_units(orbit.omega, "rad"))
          Mₛ: $(stringify_units(orbit.M_s, "M⊙"))
          aₛ: $(stringify_units(orbit.a_s, "R⊙"))
          Mₚ: $(stringify_units(orbit.M_p, "M⊙"))
          aₚ: $(stringify_units(orbit.a_p, "R⊙"))
      """
    )
end
