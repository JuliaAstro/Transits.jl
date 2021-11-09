@concrete struct SimpleOrbit <: AbstractOrbit
    period
    t0
    b
    duration
    speed
    half_period
    ref_time
end

"""
    SimpleOrbit(; period, duration, t0=0, b=0.0)

Circular orbit parameterized by the basic observables of a transiting system.

# Parameters
* `period` - The orbital period of the planets, nominally in days
* `duration` - The duration of the transit, similar units as `period`.
* `t0` - The midpoint time of the reference transit, similar units as `period`
* `b` - The impact parameter of the orbit, unitless
"""
function SimpleOrbit(;period, duration, t0=zero(period), b=0)
    half_period = 0.5 * period
    duration > half_period && error("duration cannot be longer than half the period")
    speed = 2.0 * sqrt(1.0 - b^2) / duration
    ref_time =  t0 - half_period
    SimpleOrbit(period, t0, b, duration, speed, half_period, ref_time)
end

period(orbit::SimpleOrbit) = orbit.period
duration(orbit::SimpleOrbit) = orbit.duration

relative_time(orbit::SimpleOrbit, t) = mod(t - orbit.ref_time, period(orbit)) - orbit.half_period

function relative_position(orbit::SimpleOrbit, t)
    Δt = relative_time(orbit, t)
    x = orbit.speed * Δt
    y = orbit.b
    z = abs(Δt) < 0.5 * duration(orbit) ? one(x) : -one(x)
    return SA[x, y, z]
end

# TODO: if texp, ϵ += 0.5 * texp
#function in_transit(orbit::SimpleOrbit, t)
#    Δt = relative_time.(orbit, t)
#    ϵ = 0.5 * duration(orbit)
#    return findall(x -> abs(x) < ϵ, Δt)
#end
#function in_transit(orbit::SimpleOrbit, t, r)
#    Δt = relative_time.(orbit, t)
#    ϵ = √((1.0 + r)^2.0 - orbit.b^2.0) / orbit.speed
#    return findall(x -> abs(x) < ϵ, Δt)
#end

function flip(orbit::SimpleOrbit, ror)
    t0 = orbit.t0 + orbit.half_period
    b = orbit.b / ror
    speed = orbit.speed / ror
    ref_time = orbit.t0
    return SimpleOrbit(period(orbit), t0, b, duration(orbit), speed, orbit.half_period, ref_time)
end

function Base.show(io::IO, orbit::SimpleOrbit)
    T = duration(orbit)
    P = period(orbit)
    b = orbit.b
    t0 = orbit.t0
    print(io, "SimpleOrbit(P=$P, T=$T, t0=$t0, b=$b)")
end

function Base.show(io::IO, ::MIME"text/plain", orbit::SimpleOrbit)
    T = duration(orbit)
    P = period(orbit)
    b = orbit.b
    t0 = orbit.t0
    print(io,
        "SimpleOrbit\n period: ", P,
        "\n duration: ", T,
        "\n t0: ", t0,
        "\n b: ", b
    )
end
