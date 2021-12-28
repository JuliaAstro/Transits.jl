using AstroLib: trueanom, kepler_solver
using Unitful, UnitfulAstro
using Unitful: AbstractQuantity
using KeywordCalls
using Rotations

# Domain specific unit conversions / constants/ fallbacks
Unitful.preferunits(u"Msun,Rsun,d"...)
const G_unit = Unitful.G
const G_nom = ustrip(u"Rsun^3/Msun/d^2", G_unit)

# Helpers
include("kepler_helpers.jl")

"""
    KeplerianOrbit(; kwargs...)

Keplerian orbit parameterized by the basic observables of a transiting 2-body system.

# Parameters
* `period`/`P` -- The orbital period of the planet [d].
* `t0`/`t_0` -- The midpoint time of the reference transit [d].
* `tp`/`t_p` -- The time of periastron [d].
* `duration`/`τ`/`T` -- The transit duration [d].
* `a` -- The semi-major axis [R⊙].
* `aR_star`/`aRs` -- The ratio of the semi-major axis to star radius.
* `R_planet`/`Rp` -- The radius of the planet [R⊙].
* `R_star`/`Rs` -- The radius of the star [R⊙].
* `rho_star`/`ρ_star` -- The spherical star density [M⊙/R⊙³].
* `r`/`RpRs` -- The ratio of the planet radius to star radius.
* `b` -- The impact parameter, bounded between 0 ≤ b ≤ 1.
* `ecc`/`e` -- The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1.
* `incl` -- The inclination of the orbital plane relative to the axis perpendicular to the
           reference plane [rad]
* `omega`/`ω` -- The argument of periapsis [rad].
* `cos_omega`/`cos_ω` -- The cosine of the argument of periapsis.
* `sin_omega`/`sin_ω` -- The sine of the argument of periapsis.
* `Omega`/`Ω` -- The longitude of the ascending node [rad].
* `M_planet`/`Mp` -- The mass of the planet [M⊙].
* `M_star`/`Ms` -- The mass of the star [M⊙].

# Valid combinations
The following flowchart can be used to determine which parameters can define a `KeplerianOrbit`:
1. The `period` or `a` must be given. If both given, then neither `M_star` or `rho_star` can be defined because the stellar density is now implied.
2. Only `incl` or `b` can be given.
3. If `ecc` is given, then `omega` must also be given.
4. If no stellar parameters are given, the central body is assumed to be the Sun. If only `rho_star` is given, then `R_star` is defined to be 1 solar radius. Otherwise, at most two of `M_star`, `R_star`, and `rho_star` can be given.
5. Either `t0` or `tp` must be given, but not both.
"""
@concrete struct KeplerianOrbit <: AbstractOrbit
    period
    t0
    tp
    t_ref
    duration
    a
    a_planet
    a_star
    R_planet
    R_star
    rho_planet
    rho_star
    r
    aR_star
    b
    ecc
    M0
    cos_incl
    sin_incl
    cos_omega
    sin_omega
    cos_Omega
    sin_Omega
    incl
    omega
    Omega
    n
    M_planet
    M_star
end

function KeplerianOrbit(nt::NamedTuple{(
        :period, :t0, :tp, :duration,
        :a, :R_planet, :R_star,
        :rho_star,
        :r, :aR_star, :b, :ecc, :cos_omega, :sin_omega,
        :incl, :omega, :Omega,
        :M_planet, :M_star,
    )})
    if nt.period isa Real || nt.a isa Real
        G = G_nom
    else
        G = G_unit
    end

    if (isnothing(nt.ecc) || iszero(nt.ecc)) && !isnothing(nt.duration)
        isnothing(nt.b) && throw(ArgumentError(
            "`b` must also be provided for a circular orbit if `duration given`"
        ))
        isnothing(nt.r) && throw(ArgumentError(
            "`r` must also be provided if `duration` given"
        ))
    end

    a, aR_star, period, rho_star, R_star, M_star, M_planet, duration = compute_consistent_inputs(
        nt.a, nt.aR_star, nt.period, nt.rho_star, nt.R_star, nt.M_star, nt.M_planet, G,
        nt.ecc, nt.duration, nt.b, nt.r,
    )
    r = isnothing(nt.r) ? zero(aR_star) : nt.r
    M_tot = compute_M_tot(M_star, M_planet)
    R_planet = compute_R_planet(R_star, r, nt.R_planet)
    rho_planet = compute_rho(M_planet, R_planet)

    n = 2.0 * π / period
    a_star = compute_a_X(a, M_planet, M_tot)
    a_planet = -compute_a_X(a, M_star, M_tot)

    # Omega
    Omega = nt.Omega
    sin_Omega, cos_Omega = sincos(Omega)

    # Eccentricity
    ecc = nt.ecc
    if iszero(ecc)
        M0 = 0.5 * π
        incl_factor_inv = 1.0
        omega = zero(Omega)
        cos_omega, sin_omega = 1.0, 0.0
    else
        if !isnothing(nt.omega)
            all(!isnothing, (nt.cos_omega, nt.sin_omega)) && throw(ArgumentError(
                "Only `ω`, or `cos_ω` and `sin_ω` can be provided"
            ))
            omega = nt.omega
            sin_omega, cos_omega = sincos(nt.omega)
        else
            cos_omega, sin_omega = nt.cos_omega, nt.sin_omega
            omega = atan(sin_omega, cos_omega)
        end
        E0 = compute_E0(ecc, cos_omega, sin_omega)
        M0 = compute_M0(ecc, E0)

        incl_factor_inv = compute_incl_factor_inv(ecc, sin_omega)
    end

    # Jacobian for cos(i) -> b
    dcosi_db = compute_dcosi_db(a, R_star, incl_factor_inv)

    if !isnothing(nt.b)
        if any(!isnothing, (nt.incl, duration))
            throw(ArgumentError("Only `incl`, `b`, or `duration` can be given"))
        end
        b = nt.b
        cos_incl = dcosi_db * b
        incl = acos(cos_incl)
        sin_incl = sin(incl)
        duration = nt.duration
    elseif !isnothing(nt.incl)
        !isnothing(nt.duration) && throw(ArgumentError(
            "Only `incl`, `b`, or `duration` can be given"
        ))
        incl = nt.incl
        sin_incl, cos_incl = sincos(incl)
        b = compute_b(cos_incl, dcosi_db)
        duration = nt.duration
    elseif !isnothing(nt.duration)
        duration = nt.duration
        b = compute_b(a_planet, R_star, duration, period, incl_factor_inv, ecc, sin_omega)
        cos_incl = dcosi_db * b
        incl = acos(cos_incl)
        sin_incl = sin(incl)
    else
        incl = 0.5 * π
        cos_incl = 0.0
        sin_incl = 1.0
        b = 0.0
        duration = nt.duration
    end

    # Compute remaining system parameters
    t0, tp = compute_t0_tp(nt.t0, nt.tp; M0=M0, n=n)
    t_ref = tp - t0

    # Sanitize dimensionless units
    if period isa AbstractQuantity
        r, aR_star, b = NoUnits.((r, aR_star, b))
    end

    return KeplerianOrbit(
        period, t0, tp, t_ref, duration,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        r, aR_star, b, ecc, M0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
        incl, omega, Omega,
        n,
        M_planet, M_star,
    )
end

@kwcall KeplerianOrbit(
    period=nothing, t0=nothing, tp=nothing, duration=nothing,
    a=nothing, R_planet=nothing, R_star=nothing,
    rho_star=nothing,
    r=nothing, aR_star=nothing, b=nothing, ecc=0.0, cos_omega=nothing, sin_omega=nothing,
    incl=nothing, omega=0.0, Omega=0.0,
    M_planet=nothing, M_star=nothing,
)

@kwalias KeplerianOrbit [
    P => period,
    t_0 => t0,
    t_p => tp,
    τ => duration,
    T => duration,
    aRs => aR_star,
    Rp => R_planet,
    Rs => R_star,
    ρ_star => rho_star,
    RpRs => r,
    e => ecc,
    cos_ω => cos_omega,
    sin_ω => sin_omega,
    ω => omega,
    Ω => Omega,
    Mp => M_planet,
    Ms => M_star,
]

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
    M = compute_M(t, orbit.t0, orbit.t_ref, orbit.n)
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
            tp = orbit.tp + 0.5*orbit.period,
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
            tp = orbit.tp,
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


function compute_consistent_inputs(a, aR_star, period, rho_star, R_star, M_star, M_planet, G, ecc, duration, b, r)
    if isnothing(a) && isnothing(period)
        throw(ArgumentError("At least `a` or `P` must be specified"))
    end

    if (isnothing(ecc) || iszero(ecc)) && !isnothing(duration)
        isnothing(R_star) && (R_star = compute_R_star_nom(G))
        aR_star = compute_aor(duration, period, b, r=r)
        a = compute_a(aR_star, R_star)
        duration = nothing
    end

    if isnothing(M_planet) && (!isnothing(a) || !isnothing(period))
        M_planet = compute_M_planet_nom(G)
    end

    # Compute implied stellar density
    implied_rho_star = false
    if all(!isnothing, (a, period))
        if any(!isnothing, (rho_star, M_star))
            throw(ArgumentError(
                "If both `a` and `P` are given, `rho_star` or `M_star` cannot be defined"
            ))
        end

        # Default to R_star = 1 R⊙ if not provided
        isnothing(R_star) && (R_star = oneunit(a))

        # Compute implied mass
        M_tot = compute_M_tot(a, G, period)

        # Compute implied density
        M_star = M_tot - M_planet
        rho_star = compute_rho_star(M_star, R_star)
        implied_rho_star = true
    end

    # Check combination of stellar params are valid
    if all(isnothing, (R_star, M_star))
        R_star = compute_R_star_nom(G)
        isnothing(rho_star) && (M_star = oneunit(M_planet))
    end

    if !implied_rho_star && sum(isnothing, (rho_star, R_star, M_star)) ≠ 1
        throw(ArgumentError(
            "Must provide exactly two of: `rho_star`, `R_star`, or `M_star` if rho_star not implied"
        ))
    end

    # Compute stellar params
    if isnothing(rho_star)
        rho_star = compute_rho_star(M_star, R_star)
    elseif isnothing(R_star)
        R_star = compute_R_star(rho_star, M_star)
    else
        M_star = compute_M_star(rho_star, R_star)
    end

    # Compute planet params
    M_tot = compute_M_tot(M_star, M_planet)
    if isnothing(a)
        if isnothing(aR_star)
            a = compute_a(M_tot, period, G)
            aR_star = a / R_star
        else
            a = compute_a(aR_star, R_star)
        end
    else
        period = compute_period(M_tot, a, G)
        aR_star = a / R_star
    end

    return a, aR_star, period, rho_star, R_star, M_star, M_planet, duration
end

stringify_units(value, unit_str) = @sprintf "%.4f %s" value unit_str
function stringify_units(value::Unitful.AbstractQuantity, unit_str)
    u = upreferred(value)
    return stringify_units(ustrip(u), string(unit(u)))
end
stringify_units(value::Nothing, unit) = "$(value)"
stringify_units(value) = @sprintf "%.4f"  value
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    print(
        IOContext(io, :fancy_exponent => false),
        """
        Keplerian Orbit
         P:      $(stringify_units(orbit.period, "d"))
         t₀:     $(stringify_units(orbit.t0, "d"))
         tₚ:     $(stringify_units(orbit.tp, "d"))
         t_ref:  $(stringify_units(orbit.t_ref, "d"))
         τ:      $(stringify_units(orbit.duration, "d"))
         a:      $(stringify_units(orbit.a, "R⊙"))
         aₚ:     $(stringify_units(orbit.a_planet, "R⊙"))
         aₛ:     $(stringify_units(orbit.a_star, "R⊙"))
         Rₚ:     $(stringify_units(orbit.R_planet, "R⊙"))
         Rₛ:     $(stringify_units(orbit.R_star, "R⊙"))
         ρₚ:     $(stringify_units(orbit.rho_planet, "M⊙/R⊙³"))
         ρₛ:     $(stringify_units(orbit.rho_star, "M⊙/R⊙³"))
         r:      $(stringify_units(orbit.r))
         aRₛ:    $(stringify_units(orbit.aR_star))
         b:      $(stringify_units(orbit.b))
         ecc:    $(stringify_units(orbit.ecc))
         cos(i): $(stringify_units(orbit.cos_incl))
         sin(i): $(stringify_units(orbit.sin_incl))
         cos(ω): $(stringify_units(orbit.cos_omega))
         sin(ω): $(stringify_units(orbit.sin_omega))
         cos(Ω): $(stringify_units(orbit.cos_Omega))
         sin(Ω): $(stringify_units(orbit.sin_Omega))
         i:      $(stringify_units(orbit.incl, "rad"))
         ω:      $(stringify_units(orbit.omega, "rad"))
         Ω:      $(stringify_units(orbit.Omega, "rad"))
         Mₚ:     $(stringify_units(orbit.M_planet, "M⊙"))
         Mₛ:     $(stringify_units(orbit.M_star, "M⊙"))"""
    )
end
