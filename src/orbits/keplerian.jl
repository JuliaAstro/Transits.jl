using AstroLib: trueanom, kepler_solver
using PhysicalConstants
using Unitful, UnitfulAstro
using KeywordCalls
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
    period::T
    t_0::T
    t_p::T
    t_ref::T
    a::L
    a_p::L
    a_s::L
    R_p::L
    R_s::L
    rho_p::D
    rho_s::D
    aR_s::R
    b::R
    ecc::R
    M_0::R
    incl::A
    omega::A
    Omega::A
    n::I
    M_p::M
    M_s::M
end

function normalize_inputs(
    period, t_0, t_p, t_ref,
    a, a_p, a_s, R_p, R_s,
    rho_p, rho_s,
    aR_s, b, ecc, M_0,
    incl, omega, Omega,
    n,
    M_p, M_s,
    )

    # Normalize unitless types
    aR_s, b, ecc = promote(aR_s, b, ecc)

    # Normalize quantities
    if !(period isa Real)
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
        a, a_p, a_s, R_p, R_s = uconvert.(u"Rsun", (a, a_p, a_s, R_p, R_s))
        rho_p, rho_s = uconvert.(u"Msun/Rsun^3", (rho_p, rho_s))
        incl, omega, Omega = uconvert.(u"rad", (incl, omega, Omega))
        M_p, M_s = uconvert.(u"Msun", (M_p, M_s))
    else
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
        a, a_p, a_s, R_p, R_s = promote(a, a_p, a_s, R_p, R_s)
    end

    return (
        period, t_0, t_p, t_ref,
        a, a_p, a_s, R_p, R_s,
        rho_p, rho_s,
        aR_s, b, ecc, M_0,
        incl, omega, Omega,
        n,
        M_p, M_s,
    )
end

function KeplerianOrbit(nt::NamedTuple{(
        :period, :t_0,
        :R_p, :R_s,
        :rho_s,
        :aR_s, :b, :ecc,
        :incl, :omega, :Omega,
        :M_p, :M_s,
    )})
    G = nt.period isa Real ? G_nom : G_unit
    period, t_0 = nt.period, nt.t_0
    omega, Omega = nt.omega, nt.Omega

    if !isnothing(nt.rho_s) & isnothing(nt.aR_s)
        rho_s = nt.rho_s isa Real ? convert_rho_s(nt.rho_s) : nt.rho_s
        R_s = nt.R_s
        R_p = isnothing(nt.R_p) ? zero(R_s) : nt.R_p
        aR_s, ecc = compute_aR_s(rho_s, period, G), nt.ecc
        a = compute_a(rho_s, period, R_s, G)
        M_s, a_s, M_p, a_p = compute_RV_params(rho_s, R_s, a, period, G; M_p=nt.M_p)
        rho_p = isnothing(nt.M_p) ? zero(rho_s) : compute_rho(M_p, R_p)
        if !isnothing(nt.incl) & isnothing(nt.b)
            b = compute_b(rho_s, nt.period, G, sincos(nt.incl), nt.ecc, nt.omega)
            incl = nt.incl
        elseif isnothing(nt.incl) & !isnothing(nt.b)
            b = nt.b
            incl = compute_incl(nt.rho_s, nt.period, G, b, nt.ecc, sincos(nt.omega))
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end

        # Compute remaining system parameters
        n = 2.0 * π / nt.period
        M_0 = compute_M_0(ecc, omega)
        t_p = nt.t_0 - M_0 / n
        t_ref = t_p - nt.t_0

        # Normalize inputs
        (
            period, t_0, t_p, t_ref,
            a, a_p, a_s, R_p, R_s,
            rho_p, rho_s,
            aR_s, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_p, M_s,
        ) = normalize_inputs(
            period, t_0, t_p, t_ref,
            a, a_p, a_s, R_p, R_s,
            rho_p, rho_s,
            aR_s, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_p, M_s,
        )
    elseif isnothing(nt.rho_s) & !isnothing(nt.aR_s)
        aR_s, ecc = nt.aR_s, nt.ecc
        rho_s = compute_rho_s(aR_s, period, G)
        R_s = isnothing(nt.R_s) ? (nt.period isa Real ? 1.0 : 1.0u"Rsun") : nt.R_s
        R_p = isnothing(nt.R_p) ? zero(R_s) : nt.R_p
        a = compute_a(rho_s, period, R_s, G)
        M_s, a_s, M_p, a_p = compute_RV_params(rho_s, R_s, a, period, G; M_p=nt.M_p)
        rho_p = isnothing(nt.M_p) ? zero(rho_s) : compute_rho(M_p, R_p)
        if !isnothing(nt.incl) & isnothing(nt.b)
            b = compute_b(nt.aR_s, sincos(nt.incl), nt.ecc, nt.omega)
            incl = nt.incl
        elseif isnothing(nt.incl) & !isnothing(nt.b)
            b = nt.b
            incl = compute_incl(rho_s, nt.period, G, b, nt.ecc, sincos(nt.ω))
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end

        # Compute remaining system parameters
        n = 2.0 * π / nt.period
        M_0 = compute_M_0(nt.ecc, nt.omega)
        t_p = nt.t_0 - M_0 / n
        t_ref = t_p - nt.t_0

        # Normalize inputs
        (
            period, t_0, t_p, t_ref,
            a, a_p, a_s, R_p, R_s,
            rho_p, rho_s,
            aR_s, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_p, M_s,
        ) = normalize_inputs(
            period, t_0, t_p, t_ref,
            a, a_p, a_s, R_p, R_s,
            rho_p, rho_s,
            aR_s, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_p, M_s,
        )
    else
        throw(ArgumentError("Either rho_s or aR_s must be specified"))
    end

    return KeplerianOrbit(
        period, t_0, t_p, t_ref,
        a, a_p, a_s, R_p, R_s,
        rho_p, rho_s,
        aR_s, b, ecc, M_0,
        incl, omega, Omega,
        n,
        M_p, M_s,
    )
end

@kwcall KeplerianOrbit(
    period, t_0,
    R_p=nothing, R_s=nothing,
    rho_s=nothing,
    aR_s=nothing, b=nothing, ecc,
    incl=nothing, omega, Omega,
    M_p=nothing, M_s=nothing,
)

@kwalias KeplerianOrbit [
    ρₛ => rho_s,
    aRₛ => aR_s,
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
compute_rho_s(aR_s, period, G_nom) = (3.0 * π / (G_nom * period^2)) * aR_s^3
compute_rho_s(aR_s, period, G::typeof(G_unit)) = (3.0 * π / (G * period^2)) * aR_s^3
compute_rho_s(a, period, R_s, G) = compute_rho_s(compute_aR_s(a, R_s), period, G)

# General density
compute_rho(M, R) = 0.75 * M / (π*R^3)

# Semi-major axis / star radius ratio
compute_aR_s(a, R_s) = a / R_s
compute_aR_s(rho_s, period, G_nom) = cbrt(G_nom * period^2 * rho_s / (3.0 * π))
compute_aR_s(rho_s, period, G::typeof(G_unit)) = cbrt(G * period^2 * rho_s / (3.0 * π))

# Semi-major axis
compute_a(aR_s, R_s) = aR_s * R_s
#compute_a(M_tot, period, G) = cbrt(G * M_tot * period^2 / (4.0 * π^2))
compute_a(rho_s, period, R_s, G) = compute_a(compute_aR_s(rho_s, period, G), R_s)

# Impact parameter
function compute_b(aR_s, sincos_incl, ecc, omega)
    sin_omega, cos_omega = sincos(omega)
    incl_factor_inv  = (1.0 - ecc^2) / (1.0 + ecc * sin_omega)
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
compute_M_tot(a, period, G) = 4.0 * π^2 * a^3 / (G * period^2)

function compute_RV_params(rho_s, R_s, a, period, G; M_p=nothing)
    M_tot = compute_M_tot(a, period, G)
    isnothing(M_p) && (M_p = zero(M_tot))
    M_s = M_tot - M_p
    a_p = -(M_s / M_tot) * a
    a_s = a + a_p
    return M_s, a_s, M_p, a_p
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

function compute_M_0(ecc, omega)
    sin_omega, cos_omega = sincos(omega)
    E₀ = 2.0 * atan(√(1.0 - ecc) * cos_omega, √(1.0 + ecc) * (1.0 + sin_omega))
    M_0 = E₀ - ecc * sin(E₀)
    return M_0
end

#=
function flip(orbit::KeplerianOrbit, ror)
    if iszero(orbit.ecc)
    return KeplerianOrbit(
        rho_s = orbit.M_p / ((4.0/3.0) * π * ),
        R_s = ror * orbit.R_s,
        period = orbit.period,
        ecc = orbit.ecc,
        t_0 = orbit.t_0,
        b = orbit.b,
        Omega = orbit.Omega,
        omega = orbit.omega,
        M_p = orbit.rho_s * (4.0/3.0) * π * orbit.R_s^3,
    )
    else
    return KeplerianOrbit(
        rho_s = orbit.rho_s,
        R_s = ror * orbit.R_s,
        period = orbit.period,
        ecc = orbit.ecc,
        t_0 = orbit.t_0,
        b = orbit.b,
        Omega = orbit.Omega,
        omega = orbit.omega,
        M_p = orbit.rho_s * (4.0/3.0) * π * orbit.R_s^3,
    )
    end
end
=#

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
