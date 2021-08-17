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
* `period`/`P` -- The orbital period of the planet [d].
* `t_0`/`t₀` -- The midpoint time of the reference transit [d].
* `t_p`/`tₚ` -- The time of periastron [d].
* `duration`/`τ` -- The transit duration [d].
* `a` -- The semi-major axis [R⊙].
* `R_planet`/`Rₚ` -- The radius of the planet [R⊙].
* `R_star`/`Rₛ` -- The radius of the star [R⊙].
* `rho_star`/`ρₛ` -- The spherical star density [M⊙/R⊙³].
* `RpRs` -- The ratio of the planet radius to star radius.
* `aR_star`/`aRₛ` -- The ratio of the semi-major axis to star radius.
* `b` -- The impact parameter, bounded between 0 ≤ b ≤ 1.
* `ecc` -- The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1.
* `cos_omega`/`cos_ω` -- The cosine of the argument of periapsis.
* `sin_omega`/`sin_ω` -- The sine of the argument of periapsis.
* `incl` -- The inclination of the orbital plane relative to the axis perpendicular to the
           reference plane [rad]
* `omega`/`ω` -- The argument of periapsis [rad].
* `Omega`/`Ω` -- The longitude of the ascending node [rad].
* `M_planet`/`Mₚ` -- The mass of the planet [M⊙].
* `M_star`/`Mₛ` -- The mass of the star [M⊙].

# Valid combinations
The following flowchart can be used to determine which parameters can define a `KeplerianOrbit`:
1. The `period` or `a` must be given. If both given, then neither `M_star` or `rho_star` can be defined because the stellar density is now implied.
2. Only `incl` or `b` can be given.
3. If `ecc` is given, then `omega` must also be given.
4. If no stellar parameters are given, the central body is assumed to be the Sun. If only `rho_star` is given, then `R_star` is defined to be 1 solar radius. Otherwise, at most two of `M_star`, `R_star`, and `rho_star` can be given.
5. Either `t_0` or `t_p` must be given, but not both.
"""
struct KeplerianOrbit{T,L,D,R,A,I,M} <: AbstractOrbit
    period::T
    t_0::T
    t_p::T
    t_ref::T
    duration::Union{Nothing, T}
    a::L
    a_planet::L
    a_star::L
    R_planet::Union{Nothing, L}
    R_star::L
    rho_planet::Union{Nothing, D}
    rho_star::D
    RpRs::R
    aR_star::R
    b::R
    ecc::R
    M_0::R
    cos_incl::R
    sin_incl::R
    cos_omega::R
    sin_omega::R
    cos_Omega::R
    sin_Omega::R
    incl::A
    omega::A
    Omega::A
    n::I
    M_planet::M
    M_star::M
end

function normalize_inputs(
    period, t_0, t_p, t_ref, duration,
    a, a_planet, a_star, R_planet, R_star,
    rho_planet, rho_star,
    RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
    incl, omega, Omega,
    M_planet, M_star,
    no_units,
    )

    # Normalize dimensionless quantities
    RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega = promote(
        RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega
    )

    # Normalize remaining quantities
    if no_units
        period, t_0, t_p, t_ref = promote(period, t_0, t_p, t_ref)
        !isnothing(duration) && (duration = convert(typeof(period), duration))
        a, a_planet, a_star, R_star = promote(a, a_planet, a_star, R_star)
        !isnothing(R_planet) && (R_planet = convert(typeof(a), R_planet))
        !isnothing(rho_planet) && (rho_planet = convert(typeof(rho_star), rho_planet))
        incl, omega, Omega = promote(incl, omega, Omega)
        M_planet, M_star = promote(M_planet, M_star)
    else
        period, t_0, t_p, t_ref = uconvert.(u"d", (period, t_0, t_p, t_ref))
        !isnothing(duration) && (duration = uconvert(u"d", duration))
        a, a_planet, a_star, R_star = uconvert.(u"Rsun", (a, a_planet, a_star, R_star))
        !isnothing(R_planet) && (R_planet = uconvert(u"Rsun", R_planet))
        rho_star = uconvert(u"Msun/Rsun^3", rho_star)
        !isnothing(rho_planet) && (rho_planet = uconvert(u"Msun/Rsun^3", rho_planet))
        incl, omega, Omega = uconvert.(u"rad", (incl, omega, Omega))
        M_planet, M_star = uconvert.(u"Msun", (M_planet, M_star))
    end

    return (
        period, t_0, t_p, t_ref, duration,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
        incl, omega, Omega,
        M_planet, M_star,
    )
end

function KeplerianOrbit(nt::NamedTuple{(
        :period, :t_0, :t_p, :duration,
        :a, :R_planet, :R_star,
        :rho_star,
        :RpRs, :aR_star, :b, :ecc, :cos_omega, :sin_omega,
        :incl, :omega, :Omega,
        :M_planet, :M_star,
    )})
    if nt.period isa Real || nt.a isa Real
        G = G_nom
        no_units = true
    else
        G = G_unit
        no_units = false
    end

    if (isnothing(nt.ecc) || iszero(nt.ecc)) && !isnothing(nt.duration)
        #isnothing(nt.R_star) && (R_star = no_units ? 1.0 : 1.0u"Rsun")
        isnothing(nt.b) && throw(ArgumentError(
            "`b` must also be provided for a circular orbit if `duration given`"
        ))
        isnothing(nt.RpRs) && throw(ArgumentError(
            "`RₚRₛ` must also be provided if `duration` given"
        ))
    end

    a, aR_star, period, rho_star, R_star, M_star, M_planet, duration = compute_consistent_inputs(
        nt.a, nt.aR_star, nt.period, nt.rho_star, nt.R_star, nt.M_star, nt.M_planet, G,
        nt.ecc, nt.duration, nt.b, nt.RpRs,
    )
    RpRs = isnothing(nt.RpRs) ? zero(aR_star) : nt.RpRs
    M_tot = M_star + M_planet
    R_planet = compute_R_planet(R_star, RpRs, nt.R_planet)
    rho_planet = compute_rho(M_planet, R_planet)

    n = 2.0 * π / period
    a_star = a * M_planet / M_tot
    a_planet = -a * M_star / M_tot

    # Omega
    if isnothing(nt.Omega)
        Omega = 0.0
    else
        Omega = nt.Omega
    end
    sin_Omega, cos_Omega = sincos(Omega)

    # Eccentricity
    if isnothing(nt.ecc)
        ecc = 0.0
        M_0 = 0.5 * π # TODO: find out why this fails
        incl_factor_inv = 1.0
        omega = zero(Omega)
        cos_omega, sin_omega = zero(omega), oneunit(omega)
    else
        ecc = nt.ecc
        if !isnothing(nt.omega)
            all( (!isnothing).((nt.cos_omega, nt.sin_omega)) ) && throw(ArgumentError(
                "Only `ω`, or `cos_ω` and `sin_ω` can be provided"
            ))
            omega = nt.omega
            sin_omega, cos_omega = sincos(nt.omega)
        elseif all( (!isnothing).((nt.cos_omega, nt.sin_omega)) )
            cos_omega, sin_omega = nt.cos_omega, nt.sin_omega
            omega = atan(sin_omega, cos_omega)
        else
            throw(ArgumentError("`ω` must also be provided if `ecc` specified"))
        end
        E_0 = 2.0 * atan(√(1.0 - ecc) * cos_omega, √(1.0 + ecc) * (1.0 + sin_omega))
        M_0 = E_0 - ecc * sin(E_0)

        incl_factor_inv  = (1.0 - ecc^2) / (1.0 + ecc * sin_omega)
    end

    # Jacobian for cos(i) -> b
    dcosi_db = R_star / a * (1.0 / incl_factor_inv)

    if !isnothing(nt.b)
        any( (!isnothing).((nt.incl, duration)) ) && throw(ArgumentError(
            "Only `incl`, `b`, or `duration` can be given"
        ))
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
        b = cos_incl / dcosi_db
        duration = nt.duration
    elseif !isnothing(nt.duration)
        duration = nt.duration
        c = sin(π * duration / (incl_factor_inv) / period)
        c_sq = c^2
        ecc_sin_omega = ecc*sin_omega
        aor = a_planet / R_star
        b = √(
            (aor^2 * c_sq - 1.0) /
            (
                c_sq * ecc_sin_omega^2 +
                2.0*c_sq*ecc_sin_omega +
                c_sq - ecc^4 + 2.0*ecc^2 - 1.0
            )
        ) * (1.0 - ecc^2)
        cos_incl = dcosi_db * b
        incl = acos(cos_incl)
        sin_incl = sin(incl)
    else
        incl = 0.5*π
        cos_incl = 0.0
        sin_incl = 1.0
        b = 0.0
        duration = nt.duration
    end

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

    # Normalize inputs
    (
        period, t_0, t_p, t_ref, duration,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
        incl, omega, Omega,
        M_planet, M_star,
    ) = normalize_inputs(
        period, t_0, t_p, t_ref, duration,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
        incl, omega, Omega,
        M_planet, M_star,
        no_units,
    )
    return KeplerianOrbit(
        period, t_0, t_p, t_ref, duration,
        a, a_planet, a_star, R_planet, R_star,
        rho_planet, rho_star,
        RpRs, aR_star, b, ecc, M_0, cos_incl, sin_incl, cos_omega, sin_omega, cos_Omega, sin_Omega,
        incl, omega, Omega,
        n,
        M_planet, M_star,
    )
end

@kwcall KeplerianOrbit(
    period=nothing, t_0=nothing, t_p=nothing, duration=nothing,
    a=nothing, R_planet=nothing, R_star=nothing,
    rho_star=nothing,
    RpRs=nothing, aR_star=nothing, b=nothing, ecc=nothing, cos_omega=nothing, sin_omega=nothing,
    incl=nothing, omega=nothing, Omega=nothing,
    M_planet=nothing, M_star=nothing,
)

@kwalias KeplerianOrbit [
    P => period,
    t₀ => t_0,
    tₚ => t_p,
    τ => duration,
    Rₚ => R_planet,
    Rₛ => R_star,
    ρₛ => rho_star,
    aRₛ => aR_star,
    cos_ω => cos_omega,
    sin_ω => sin_omega,
    ω => omega,
    Ω => Omega,
    Mₚ => M_planet,
    Mₛ => M_star,
]

# `M_planet` << `M_star` approxs
#compute_rho_star(aR_star, period, G_nom) = (3.0 * π / (G_nom * period^2)) * aR_star^3
#compute_rho_star(aR_star, period, G::typeof(G_unit)) = (3.0 * π / (G * period^2)) * aR_star^3
#compute_rho_star(a, period, R_star, G) = compute_rho_star(compute_aR_star(a, R_star), period, G)
#compute_aR_star(a, R_star) = a / R_star
#compute_aR_star(rho_star, period, G_nom) = cbrt(G_nom * period^2 * rho_star / (3.0 * π))
#compute_aR_star(rho_star, period, G::typeof(G_unit)) = cbrt(G * period^2 * rho_star / (3.0 * π))

# Spherical density
compute_rho(M, R) = 0.75 * M / (π*R^3)
compute_rho(M, R::Nothing) = nothing

# Semi-major axis / star radius ratio, assuming circular orbit
function compute_aor(duration, period, b; RpRs=nothing)
    RpRs = isnothing(RpRs) ? 0.0 : RpRs
    sin_ϕ, cos_ϕ = sincos(π * duration / period)
    return √( (1 + RpRs)^2 - (b*cos_ϕ)^2 ) / sin_ϕ
end

# Planet radius
compute_R_planet(R_star, RpRs, R_planet) = R_planet
compute_R_planet(R_star, RpRs, R_planet::Nothing) = iszero(RpRs) ? nothing : R_star * RpRs

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

function compute_consistent_inputs(
    a, aR_star, period, rho_star, R_star, M_star, M_planet,
    G, ecc, duration, b, RpRs
    )
    all( isnothing.((a, period)) ) && throw(
        ArgumentError("At least `a` or `P` must be specified")
    )

    no_units = G isa Real

    if (isnothing(ecc) || iszero(ecc)) && !isnothing(duration)
        isnothing(R_star) && (R_star = no_units ? 1.0 : 1.0u"Rsun")
        aR_star = compute_aor(duration, period, b, RpRs=RpRs)
        a = R_star * aR_star
        duration = nothing
    end

    if !isnothing(a) && isnothing(M_planet)
        M_planet = no_units ? 0.0 : 0.0u"Msun"
    end

    if !isnothing(period) && isnothing(M_planet)
        M_planet = no_units ? 0.0 : 0.0u"Msun"
    end

    # Compute implied stellar density
    implied_rho_star = false
    if all( (!isnothing).((a, period)) )
        if any( (!isnothing).((rho_star, M_star)) )
            throw(ArgumentError(
                "If both `a` and `P` are given, `ρₛ` or `Mₛ` cannot be defined"
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
            "Must provide exactly two of: `ρₛ`, `Rₛ`, or `Mₛ` if ρₛ not implied"
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

    return a, aR_star, period, rho_star, R_star, M_star, M_planet, duration
end

stringify_units(value::Unitful.AbstractQuantity, unit) = value
stringify_units(value, unit) = "$value $unit"
function Base.show(io::IO, ::MIME"text/plain", orbit::KeplerianOrbit)
    if orbit.period isa Real
        if isnothing(orbit.rho_planet)
            rho_planet = nothing
        else
            rho_planet = orbit.rho_planet * MsunRsun_to_gcc
        end
        rho_star = orbit.rho_star * MsunRsun_to_gcc
    else
        if isnothing(orbit.rho_planet)
            rho_planet = nothing
        else
            println("we here")
            rho_planet = orbit.rho_planet |> u"g/cm^3"
        end
        rho_star = orbit.rho_star |> u"g/cm^3"
    end
    print(
        io,
        """
        Keplerian Orbit
         P: $(stringify_units(orbit.period, "d"))
         t₀: $(stringify_units(orbit.t_0, "d"))
         tₚ: $(stringify_units(orbit.t_p, "d"))
         t_ref: $(stringify_units(orbit.t_ref, "d"))
         τ: $(stringify_units(orbit.duration, "d"))
         a: $(stringify_units(orbit.a, "R⊙"))
         aₚ: $(stringify_units(orbit.a_planet, "R⊙"))
         aₛ: $(stringify_units(orbit.a_star, "R⊙"))
         Rₚ: $(stringify_units(orbit.R_planet, "R⊙"))
         Rₛ: $(stringify_units(orbit.R_star, "R⊙"))
         ρₚ: $(stringify_units(rho_planet, "g/cm³"))
         ρₛ: $(stringify_units(rho_star, "g/cm³"))
         RpRs: $(orbit.RpRs)
         aRₛ: $(orbit.aR_star)
         b: $(orbit.b)
         ecc: $(orbit.ecc)
         cos(incl): $(orbit.cos_incl)
         sin(incl): $(orbit.sin_incl)
         cos(omega): $(orbit.cos_omega)
         sin(omega): $(orbit.sin_omega)
         cos(Omega): $(orbit.cos_Omega)
         sin(Omega): $(orbit.sin_Omega)
         incl: $(stringify_units(orbit.incl, "rad"))
         ω: $(stringify_units(orbit.omega, "rad"))
         Ω: $(stringify_units(orbit.Omega, "rad"))
         Mₚ: $(stringify_units(orbit.M_planet, "M⊙"))
         Mₛ: $(stringify_units(orbit.M_star, "M⊙"))"""
    )
end
