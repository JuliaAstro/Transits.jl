using AstroLib: trueanom, kepler_solver
using PhysicalConstants
using Unitful, UnitfulAstro
using KeywordCalls
using Rotations

# Domain specific unit conversions / Constants
const G_nom = 2942.2062175044193 # Rsun^3/Msun/d^2
const G_unit = PhysicalConstants.CODATA2018.G

"""
    KeplerianOrbit(; kwargs...)

Keplerian orbit parameterized by the basic observables of a transiting 2-body system.

# Parameters
* `a` - The semi-major axis, nominally in AU
* `aR_star`/`aRₛ` - The ratio of the semi-major axis to the star radius.
* `b` - The impact parameter, bounded between 0 ≤ b ≤ 1
* `ecc` - The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1
* `period`/`P` - The orbital period of the planet, nominally in days.
* `rho_star`/`ρₛ` - The spherical star density, nominally in g/cm^3.
* `R_star`/`Rₛ` - The star mass, nominally in solar radii.
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
    a_planet::L
    a_star::L
    R_planet::L
    R_star::L
    rho_planet::D
    rho_star::D
    aR_star::R
    b::R
    ecc::R
    M_0::R
    incl::A
    omega::A
    Omega::A
    n::I
    M_planet::M
    M_star::M
end

function normalize_inputs(
    period, t_0, t_p, t_ref,
    a, a_planet, a_star, R_planet, R_star,
    rho_planet, rho_star,
    aR_star, b, ecc, M_0,
    incl, omega, Omega,
    n,
    M_planet, M_star,
    )

    # Normalize unitless types
    aR_star, b, ecc = promote(aR_star, b, ecc)

    # Normalize quantities
    if !(period isa Real)
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
        a, a_planet, a_star, R_planet, R_star = uconvert.(u"Rsun", (a, a_planet, a_star, R_planet, R_star))
        rho_planet, rho_star = uconvert.(u"Msun/Rsun^3", (rho_planet, rho_star))
        incl, omega, Omega = uconvert.(u"rad", (incl, omega, Omega))
        M_planet, M_star = uconvert.(u"Msun", (M_planet, M_star))
    else
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
        a, a_planet, a_star, R_planet, R_star = promote(a, a_planet, a_star, R_planet, R_star)
    end

    return (
        period, t_0, t_p, t_ref,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        n,
        M_planet, M_star,
    )
end

function KeplerianOrbit(nt::NamedTuple{(
        :period, :t_0,
        :a, :R_planet, :R_star,
        :rho_star,
        :aR_star, :b, :ecc,
        :incl, :omega, :Omega,
        :M_planet, :M_star,
    )})
    G = nt.period isa Real ? G_nom : G_unit
    a, period, rho_star, R_star, M_star, M_planet, G = compute_consistent_inputs(
        nt.a, nt.period, nt.rho_star, nt.R_star, nt.M_star, nt.M_planet, G
    )
    M_tot = M_star + M_planet

    n = 2.0 * π / period
    a_star = a * M_planet / M_tot
    a_planet = -a * M_star / M_tot

    # Eccentricity
    if any(isnothing.((ecc, nt.omega)))
        throw(ArgumentError(both `e` and `ω` must be provided))
    else
        omega = nt.omega
        sin_omega, cos_omega = sincos(omega)
    end

    if isnothing(nt.ecc) || iszero(nt.ecc)
        ecc = 0.0
        M_0 = 0.5*π
        incl_factor_inv = 1.0
    else
        ecc = nt.ecc
        E_0 = 2.0 * atan(√(1.0 - ecc) * cos_omega, √(1.0 + ecc) * (1.0 + sin_omega))
        M_0 = E_0 - ecc * sin(E_0)
        incl_factor_inv  = (1.0 - ecc^2) / (1.0 + ecc * sin_omega)
    end

    if !isnothing(nt.incl) & isnothing(nt.b)
        b = compute_b(nt.aR_star, sincos(nt.incl), nt.ecc, nt.omega)
        incl = nt.incl
    elseif isnothing(nt.incl) & !isnothing(nt.b)
        b = nt.b
        incl = compute_incl(rho_star, nt.period, G, b, nt.ecc, sincos(nt.ω))
    else
        throw(ArgumentError("Either incl or b must be specified"))
    end




    if !isnothing(nt.rho_star) & isnothing(nt.aR_star)
        rho_star = nt.rho_star
        R_star = nt.R_star
        R_planet = isnothing(nt.R_planet) ? 0.01*oneunit(R_star) : nt.R_planet
        aR_star, ecc = compute_aR_star(rho_star, period, G), nt.ecc
        a = compute_a(rho_star, period, R_star, G)
        M_star = rho_star * (4.0/3.0)*π*R_star^3
        a_star, M_planet, a_planet = compute_RV_params(a, M_star; M_planet=nt.M_planet)
        rho_planet = isnothing(nt.M_planet) ? zero(rho_star) : compute_rho(M_planet, R_planet)
        if !isnothing(nt.incl) & isnothing(nt.b)
            b = compute_b(rho_star, nt.period, G, sincos(nt.incl), nt.ecc, nt.omega)
            incl = nt.incl
        elseif isnothing(nt.incl) & !isnothing(nt.b)
            b = nt.b
            incl = compute_incl(nt.rho_star, nt.period, G, b, nt.ecc, sincos(nt.omega))
        else
            throw(ArgumentError("Either incl or b must be specified"))
        end

        # Compute remaining system parameters
        t_p = nt.t_0 - M_0 / n
        t_ref = t_p - nt.t_0

        # Normalize inputs
        (
            period, t_0, t_p, t_ref,
            a, a_planet, a_star, R_planet, R_star,
            rho_planet, rho_star,
            aR_star, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_planet, M_star,
        ) = normalize_inputs(
            period, t_0, t_p, t_ref,
            a, a_planet, a_star, R_planet, R_star,
            rho_planet, rho_star,
            aR_star, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_planet, M_star,
        )
    elseif isnothing(nt.rho_star) & !isnothing(nt.aR_star)
        aR_star, ecc = nt.aR_star, nt.ecc
        rho_star = compute_rho_star(aR_star, period, G)
        R_star = isnothing(nt.R_star) ? (nt.period isa Real ? 1.0 : 1.0u"Rsun") : nt.R_star
        R_planet = isnothing(nt.R_planet) ? zero(R_star) : nt.R_planet
        a = compute_a(rho_star, period, R_star, G)
        M_star = rho_star * (4.0/3.0) * π * R_star^3
        a_star, M_planet, a_planet = compute_RV_params(a, M_star; M_planet=nt.M_planet)
        rho_planet = isnothing(nt.M_planet) ? zero(rho_star) : compute_rho(M_planet, R_planet)

        # Compute remaining system parameters
        n = 2.0 * π / nt.period
        M_0 = compute_M_0(nt.ecc, nt.omega)
        t_p = nt.t_0 - M_0 / n
        t_ref = t_p - nt.t_0

        # Normalize inputs
        (
            period, t_0, t_p, t_ref,
            a, a_planet, a_star, R_planet, R_star,
            rho_planet, rho_star,
            aR_star, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_planet, M_star,
        ) = normalize_inputs(
            period, t_0, t_p, t_ref,
            a, a_planet, a_star, R_planet, R_star,
            rho_planet, rho_star,
            aR_star, b, ecc, M_0,
            incl, omega, Omega,
            n,
            M_planet, M_star,
        )
    else
        throw(ArgumentError("Either rho_star or aR_star must be specified"))
    end

    return KeplerianOrbit(
        period, t_0, t_p, t_ref,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        n,
        M_planet, M_star,
    )
end

@kwcall KeplerianOrbit(
    period=nothing, t_0=nothing,
    a=nothing, R_planet=nothing, R_star=nothing,
    rho_star=nothing,
    aR_star=nothing, b=nothing, ecc=nothing,
    incl=nothing, omega=nothing, Omega=nothing,
    M_planet=nothing, M_star=nothing,
)

@kwalias KeplerianOrbit [
    ρₛ => rho_star,
    aRₛ => aR_star,
    Rₛ => R_star,
    P => period,
    t₀ => t_0,
    Ω => Omega,
    ω => omega,
]

#############
# Orbit logic
#############
# Star density
compute_rho_star(aR_star, period, G_nom) = (3.0 * π / (G_nom * period^2)) * aR_star^3
compute_rho_star(aR_star, period, G::typeof(G_unit)) = (3.0 * π / (G * period^2)) * aR_star^3
compute_rho_star(a, period, R_star, G) = compute_rho_star(compute_aR_star(a, R_star), period, G)

# General density
compute_rho(M, R) = 0.75 * M / (π*R^3)

# Semi-major axis / star radius ratio
compute_aR_star(a, R_star) = a / R_star
compute_aR_star(rho_star, period, G_nom) = cbrt(G_nom * period^2 * rho_star / (3.0 * π))
compute_aR_star(rho_star, period, G::typeof(G_unit)) = cbrt(G * period^2 * rho_star / (3.0 * π))

# Semi-major axis
compute_a(aR_star, R_star) = aR_star * R_star
compute_a(M_tot, period, G) = cbrt(G * M_tot * period^2 / (4.0 * π^2))
compute_a(rho_star, period, R_star, G) = compute_a(compute_aR_star(rho_star, period, G), R_star)

# Impact parameter
function compute_b(aR_star, sincos_incl, ecc, omega)
    sin_omega, cos_omega = sincos(omega)
    cos_incl = sincos_incl[2]
    incl_factor_inv  = (1.0 - ecc^2) / (1.0 + ecc * sin_omega)
    return aR_star * cos_incl * incl_factor_inv
end
compute_b(rho_star, period, G, sincos_incl, ecc, omega) = compute_b(compute_aR_star(rho_star, period, G), sincos_incl, ecc, omega)

# Inclination
function compute_incl(aR_star, b, ecc, sincosomega)
    return acos((b/aR_star) * (1.0 + ecc*sincosomega[1])/(1.0 - ecc^2))
end
compute_incl(rho_star, period, G, b, ecc, sincosomega) = compute_incl(compute_aR_star(rho_star, period, G), b, ecc, sincosomega)

###########
# RV params
###########
function compute_RV_params(a, M_star; M_planet=nothing)
    isnothing(M_planet) && (M_planet = zero(M_star))
    M_tot = M_star + M_planet
    a_planet = -(M_star / M_tot) * a
    a_star = a + a_planet
    return a_star, M_planet, a_planet
end

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
# separation: aR_star, a_star / R_star, or a_planet / R_star
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
_star_position(orb, R_star, t) = _position.(orb, orb.a_star / R_star, t)
_planet_position(orb, R_star, t) = _position.(orb, orb.a_planet / R_star, t)
relative_position(orbit::KeplerianOrbit, t) = _position(orbit, -orbit.aR_star, t)

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
    end
    return M_0
end

function flip(orbit::KeplerianOrbit, ror)
    R_planet = ror * orbit.R_star
    M_planet = orbit.M_planet
    if iszero(orbit.ecc)
        return KeplerianOrbit(
            rho_star = compute_rho(M_planet, R_planet),
            R_star = R_planet,
            period = orbit.period,
            ecc = orbit.ecc,
            t_0 = orbit.t_0,
            incl = orbit.incl,
            Omega = orbit.Omega,
            omega = orbit.omega,
            M_planet = orbit.M_star,
        )
    else
        return KeplerianOrbit(
            rho_star = compute_rho(M_planet, R_planet),
            R_star = R_planet,
            period = orbit.period,
            ecc = orbit.ecc,
            t_0 = orbit.t_0,
            incl = orbit.incl,
            Omega = orbit.Omega,
            omega = orbit.omega - π,
            M_planet = orbit.M_star,
            R_planet = orbit.R_star,
        )
    end
end

stringify_units(value::Unitful.AbstractQuantity, unit) = value
stringify_units(value, unit) = "$value $unit"
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    print(
        io,
        """
        Keplerian Orbit
         a: $(stringify_units(orbit.a, "R⊙"))
         aRₛ: $(orbit.aR_star)
         b: $(orbit.b)
         ecc: $(orbit.ecc)
         P: $(stringify_units(orbit.period, "d"))
         ρₚ: $(stringify_units(orbit.rho_planet, "M⊙/R⊙³"))
         ρₛ: $(stringify_units(orbit.rho_star, "M⊙/R⊙³"))
         Rₚ: $(stringify_units(orbit.R_planet, "R⊙"))
         Rₛ: $(stringify_units(orbit.R_star, "R⊙"))
         t₀: $(stringify_units(orbit.t_0, "d"))
         tₚ: $(stringify_units(orbit.t_p, "d"))
         t_ref: $(stringify_units(orbit.t_ref, "d"))
         incl: $(stringify_units(orbit.incl, "rad"))
         Ω: $(stringify_units(orbit.Omega, "rad"))
         ω: $(stringify_units(orbit.omega, "rad"))
         Mₛ: $(stringify_units(orbit.M_star, "M⊙"))
         aₛ: $(stringify_units(orbit.a_star, "R⊙"))
         Mₚ: $(stringify_units(orbit.M_planet, "M⊙"))
         aₚ: $(stringify_units(orbit.a_planet, "R⊙"))"""
    )
end

function compute_consistent_inputs(a, aR_star, period, rho_star, R_star, M_star, M_planet, G)
    all(isnothing.((a, period))) && throw(
        ArgumentError("at least `a` or `P` must be specified")
    )

    isnothing(M_planet) && (M_planet = G isa Real ? 0.0 : 0.0u"Msun")

    # Compute implied stellar density
    implied_rho_star = false
    if all((!isnothing).((a, period)))
        if any((!isnothing).((rho_star, M_star)))
            throw(ArgumentError(
                "if both `a` and `P` are given,
                `ρₛ` or `Mₛ` cannot be defined"
            ))
        end

        # Default to Rₛ = 1 R⊙ if not provided
        isnothing(R_star) && (R_star = oneunit(a))

        # Compute implied mass
        M_tot = 4.0 * π^2 * a^3 / (G * period^2)

        # Compute implied density
        M_star = M_tot - M_planet
        rho_star = M_star / ((4.0/3.0) * π * R_star^3)
        implied_rho_star = true
    end

    # Check combination of stellar params are valid
    if all(isnothing.((R_star, M_star)))
        R_star = G isa Real ? 1.0 : 1.0u"Rsun"
        isnothing(rho_star) && (M_star = oneunit(M_planet))
    end

    if !implied_rho_star && sum(isnothing.((rho_star, R_star, M_star))) ≠ 1
        throw(ArgumentError(
            "values mut be provided for exactly two of: `ρₛ`, `Rₛ`, `Mₛ`"
        ))
    end

    # Compute stellar params
    if isnothing(rho_star)
        rho_star = 3.0 * M_star / (4.0 * π * R_star^3)
    elseif isnothing(R_star)
        R_star = ( 3.0 * M_star / (4.0 * π * rho_star) )^(1/3)
    else
        M_star = 4.0 * π * R_star^3 * rho_star / 3.0
    end

    # Compute planet params
    M_tot = M_star + M_planet
    if isnothing(a)
        a = isnothing(aR_star) ? ( G * M_tot * period^2 / (4.0 * π^2) )^(1/3) : aR_star * R_star
    else
        period = 2.0 * π * a^(3/2) / (√(G * M_tot))
    end

    return a, period, rho_star, R_star, M_star, M_planet
end
