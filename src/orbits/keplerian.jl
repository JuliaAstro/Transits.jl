using AstroLib: trueanom, kepler_solver
using PhysicalConstants
using Unitful, UnitfulAstro
using KeywordCalls
using Rotations

# Domain specific unit conversions / Constants
const G_nom = 2942.2062175044193 # Rsun^3/Msun/d^2
const G_unit = PhysicalConstants.CODATA2018.G
const MsunRsun_to_gcc = (1.0u"Msun/Rsun^3" |> u"g/cm^3").val

"""
    KeplerianOrbit(; kwargs...)

Keplerian orbit parameterized by the basic observables of a transiting 2-body system.

# Parameters
* `a` - The semi-major axis, nominally in AU
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
    n::I
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
    cos_omega::A
    sin_omega::A
    Omega::A
    cos_Omega::A
    sin_Omega::A
    M_planet::M
    M_star::M
end

function normalize_inputs(
    period, t_0, t_p, t_ref,
    n,
    a, a_planet, a_star, R_planet, R_star,
    rho_planet, rho_star,
    aR_star, b, ecc, M_0,
    incl, omega, Omega,
    M_planet, M_star,
    no_units,
    )

    # Normalize unitless types
    aR_star, b, ecc, M_0 = promote(aR_star, b, ecc, M_0)

    # Normalize quantities
    if no_units
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
        a, a_planet, a_star, R_planet, R_star = promote(a, a_planet, a_star, R_planet, R_star)
    else
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
        a, a_planet, a_star, R_planet, R_star = uconvert.(u"Rsun", (a, a_planet, a_star, R_planet, R_star))
        rho_planet, rho_star = uconvert.(u"Msun/Rsun^3", (rho_planet, rho_star))
        incl, omega, Omega = uconvert.(u"rad", (incl, omega, Omega))
        M_planet, M_star = uconvert.(u"Msun", (M_planet, M_star))
    end

    return (
        period, t_0, t_p, t_ref,
        n,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        M_planet, M_star,
    )
end

function KeplerianOrbit(nt::NamedTuple{(
        :period, :t_0, :t_p,
        :a, :R_planet, :R_star,
        :rho_star,
        :aR_star, :b, :ecc,
        :incl, :omega, :cos_Omega, :sin_Omega, :Omega, :cos_Omega, :sin_Omega,
        :M_planet, :M_star,
    )})
    if nt.period isa Real || nt.a isa Real
        G = G_nom
        no_units = true
    else
        G = G_unit
        no_units = false
    end

    if isnothing(nt.ecc) && !isnothing(nt.duration)
        isnothing(nt.R_star) && (R_star = one())
        isnothing(nt.b) && throw(ArgumentError(
            "`b` must also be provided for a circular orbit if `duration given`"
        ))
        isnothing(nt.ror) %% throw(ArgumentError(
            "`ror` must also be provided if `duration` given"
        ))
        aR_star = compute_aR_star(nt.duration, nt.period, nt.b, ror=ror)
        a = nt.R_star * aR_star
    end

    a, period, rho_star, R_star, M_star, M_planet = compute_consistent_inputs(
        nt.a, nt.period, nt.rho_star, nt.R_star, nt.M_star, nt.M_planet, G,
    )
    M_tot = M_star + M_planet
    if isnothing(nt.R_planet)
        R_planet = 0.001*oneunit(R_star)
    else
        R_planet = nt.R_planet
    end
    rho_planet = compute_rho(M_planet, R_planet)

    n = 2.0 * π / period
    a_star = a * M_planet / M_tot
    a_planet = -a * M_star / M_tot

    # Omega
    if isnothing(nt.Omega)
        Omega = nothing
    else
        Omega = nt.Omega
        sin_Omega, cos_Omega = sincos(nt.Omega)
    end

    # Eccentricity
    if isnothing(nt.ecc)
        ecc = nothing
        M_0 = 0.5 * π # TODO: find out why this fails
        incl_factor_inv = 1.0
    else
        ecc = nt.ecc

        if !isnothing(nt.omega)
            all( (!isnothing).(nt.cos_omega, nt.sin_omega) ) && throw(ArgumentError(
                "Only `cos_ω` or `sin_ω` can be provided"
            ))
            omega = nt.omega
            sin_omega, cos_omega = sincos(nt.omega)
        elseif all( (!isnothing).(cos_omega, sin_omega) )
            cos_omega, sin_omega = nt.cos_omega, nt.sin_omega
            omega = atan(sin_omega, cos_omega)
        else
            throw(ArgumentError("both `e` and `ω` must be provided"))
        end

        E_0 = 2.0 * atan(√(1.0 - ecc) * cos_omega, √(1.0 + ecc) * (1.0 + sin_omega))
        M_0 = E_0 - ecc * sin(E_0)

        incl_factor_inv  = (1.0 - ecc^2) / (1.0 + ecc * sin_omega)
    end

    # Jacobian for cos(i) -> b
    dcosi_db = R_star / a * (1.0 / incl_factor_inv)
    if !isnothing(nt.b)
        any( (!isnothing.(nt.incl, nt.duration)) ) && throw(ArgumentError(
            "Only `incl`, `b`, or `duration` can be given"
        ))
        b = nt.b
        cos_incl = dcosi_db * b
        incl = acos(cos_incl)
    elseif !isnothing(nt.incl)
        !isnothing(nt.duration) && throw(ArgumentError(
            "Only `incl`, `b`, or `duration` can be given"
        ))
        incl = nt.incl
        sin_incl, cos_incl = sincos(incl)
        b = cos_incl / dcosi_db
    elseif !isnothing(nt.duration)
        !isnothing(ecc) && throw(ArgumentError(
            "Fitting with `duration` only works for eccentric orbits"
        ))
        duration = nt.duration
        c = sin(π * duration / (incl_factor_inv) / period)
        aR_star =
    else
    end
    #=
    if !isnothing(nt.incl) & isnothing(nt.b)
        incl = nt.incl
        b = compute_b(a, R_star, sincos(incl), ecc, omega, incl_factor_inv)
    elseif isnothing(nt.incl) & !isnothing(nt.b)
        b = nt.b
        incl = compute_incl(a, R_star, b, ecc, sincos(omega))
    else
    end
    =#

    # Compute remaining system parameters
    !(isnothing(nt.t_0) ⊻ isnothing(nt.t_p)) && throw(
        ArgumentError("Either t₀ or tₚ must be specified")
    )
    if isnothing(nt.t_0)
        t_p = nt.t_p
        t_0 = t_p + M_0 / n
    else
        t_0 = nt.t_0
        t_p = t_0 - M_0 / n
    end

    t_ref = t_p - t_0

    Omega = isnothing(nt.Omega) ? (no_units ? 0.5*π : 0.5*π*u"rad") : nt.Omega

    # Normalize inputs
    (
        period, t_0, t_p, t_ref,
        n,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        M_planet, M_star,
    ) = normalize_inputs(
        period, t_0, t_p, t_ref,
        n,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        M_planet, M_star,
        no_units,
    )
    return KeplerianOrbit(
        period, t_0, t_p, t_ref,
        n,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        aR_star, b, ecc, M_0,
        incl, omega, Omega,
        M_planet, M_star,
    )
end

@kwcall KeplerianOrbit(
    period=nothing, t_0=nothing, t_p=nothing,
    a=nothing, R_planet=nothing, R_star=nothing,
    rho_star=nothing,
    aR_star=nothing, b=nothing, ecc=nothing,
    incl=nothing, omega=nothing, Omega=nothing,
    M_planet=nothing, M_star=nothing,
)

@kwalias KeplerianOrbit [
    aRₛ => aR_star,
    ρₛ => rho_star,
    Rₛ => R_star,
    P => period,
    t₀ => t_0,
    tₚ => t_p,
    Ω => Omega,
    ω => omega,
]

#############
# Orbit logic
#############
# Star density
#compute_rho_star(aR_star, period, G_nom) = (3.0 * π / (G_nom * period^2)) * aR_star^3
#compute_rho_star(aR_star, period, G::typeof(G_unit)) = (3.0 * π / (G * period^2)) * aR_star^3
#compute_rho_star(a, period, R_star, G) = compute_rho_star(compute_aR_star(a, R_star), period, G)

# Spherical density
compute_rho(M, R) = 0.75 * M / (π*R^3)

# Semi-major axis / star radius ratio
function compute_aR_star(duration, period, b; ror=Nothing)
    ror = isnothing(ror) ? zero(b) : ror
    sin_ϕ, cos_ϕ = sincos(duration / period)
    return √((1 + ror)^2 - b^2*cos_ϕ^2) / sin_ϕ
end
#compute_aR_star(a, R_star) = a / R_star
#compute_aR_star(rho_star, period, G_nom) = cbrt(G_nom * period^2 * rho_star / (3.0 * π))
#compute_aR_star(rho_star, period, G::typeof(G_unit)) = cbrt(G * period^2 * rho_star / (3.0 * π))

# Semi-major axis
# compute_a(aR_star, R_star) = aR_star * R_star
compute_a(M_tot, period, G) = cbrt(G * M_tot * period^2 / (4.0 * π^2))

# Impact parameter
function compute_b(aR_star, sincos_incl, ecc, omega, incl_factor_inv)
    sin_omega, cos_omega = sincos(omega)
    cos_incl = sincos_incl[2]
    return (aR_star) * cos_incl * incl_factor_inv
end
compute_b(a, R_star, sincos_incl, ecc, omega, incl_factor_inv) = compute_b(a/R_star, sincos_incl, ecc, omega, incl_factor_inv)

# Inclination
function compute_incl(aR_star, b, ecc, sincosomega)
    return acos((b/aR_star) * (1.0 + ecc*sincosomega[1])/(1.0 - ecc^2))
end
compute_incl(a, R_star, b, ecc, sincosomega) = compute_incl(a/R_star, b, ecc, sincosomega)

# Finds the position `r` of the planet along its orbit after rotating
# through the true anomaly `ν`, then transforms this from the
# orbital plan to the equatorial plane
# separation: aR_star, a/R_star, a_star / R_star, or a_planet / R_star
# TODO: consider moving this to a separate orbital calculations package in the future
function _position(orbit, separation, t)
    sin_ν, cos_ν = compute_true_anomaly(orbit, t)
    if isnothing(orbit.ecc) || iszero(orbit.ecc)
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
    if isnothing(orbit.ecc) || iszero(orbit.ecc)
        return sincos(M)
    else
        E = kepler_solver(M, orbit.ecc)
        return sincos(trueanom(E, orbit.ecc))
    end
end

function flip(orbit::KeplerianOrbit, R_planet)
    M_planet = orbit.M_planet
    if isnothing(orbit.ecc) || iszero(orbit.ecc)
        return KeplerianOrbit(
            period = orbit.period,
            t_p = orbit.t_p + 0.5*orbit.period,
            incl = orbit.incl,
            Omega = orbit.Omega,
            omega = orbit.omega,
            M_star = orbit.M_planet,
            M_planet = orbit.M_star,
            R_star = R_planet,
            R_planet = orbit.R_star,
            ecc = orbit.ecc,
        )
    else
        return KeplerianOrbit(
            period = orbit.period,
            t_p = orbit.t_p,
            incl = orbit.incl,
            Omega = orbit.Omega,
            omega = orbit.omega - π,
            M_star = orbit.M_planet,
            M_planet = orbit.M_star,
            R_star = R_planet,
            R_planet = orbit.R_star,
            ecc = orbit.ecc,
        )
    end
end

function compute_consistent_inputs(a, period, rho_star, R_star, M_star, M_planet, G)
    all( isnothing.((a, period)) ) && throw(
        ArgumentError("at least `a` or `P` must be specified")
    )

    no_units = G isa Real
    if !isnothing(a) && isnothing(M_planet)
        M_planet = no_units ? 0.0 : 0.0u"Msun"
    end

    # Compute implied stellar density
    implied_rho_star = false
    if all( (!isnothing).((a, period)) )
        if any( (!isnothing).((rho_star, M_star)) )
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
    if all( isnothing.((R_star, M_star)) )
        R_star = no_units ? 1.0 : 1.0u"Rsun"
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
        if isnothing(aR_star)
            a = ( G * M_tot * period^2 / (4.0 * π^2) )^(1/3)
            aR_star = a / R_star
        else
            a = aR_star * R_star
        end
    else
        period = 2.0 * π * a^(3/2) / (√(G * M_tot))
        aR_star = a / R_star
    end

    return a, period, rho_star, R_star, M_star, M_planet
end

stringify_units(value::Unitful.AbstractQuantity, unit) = value
stringify_units(value, unit) = "$value $unit"
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    if orbit.period isa Real
        rho_planet = orbit.rho_planet * MsunRsun_to_gcc
        rho_star = orbit.rho_star * MsunRsun_to_gcc
    else
        rho_planet = orbit.rho_planet |> u"g/cm^3"
        rho_star = orbit.rho_star |> u"g/cm^3"
    end
    print(
        io,
        """
        Keplerian Orbit
         a: $(stringify_units(orbit.a, "R⊙"))
         aRₛ: $(orbit.aR_star)
         b: $(orbit.b)
         ecc: $(orbit.ecc)
         P: $(stringify_units(orbit.period, "d"))
         ρₚ: $(stringify_units(rho_planet, "M⊙/R⊙³"))
         ρₛ: $(stringify_units(rho_star, "g/cm³"))
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
