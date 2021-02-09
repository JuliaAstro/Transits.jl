struct SimpleOrbit <: AbstractOrbit
    period
    t0
    b
    duration
    r_star

    b_norm
    speed
    half_period
    ref_time
end

"""
    SimpleOrbit(; period, duration, t0=0, b=0, r_star=1)

Circular orbit parameterized by the basic observables of a transiting system.

# Parameters
* `period` - The orbital period of the planets, nominally in days
* `duration` The duration of the transit, same units as `period`
* `t0` - The midpoint time of the reference transit, same units as `period`
* `b` - The impact parameter of the orbit
* `r_star` - The radius of the star, nominally in solar radii
"""
function SimpleOrbit(;period, duration, t0=zero(period), b=0.0, r_star=1.0)
    half_period = 0.5 * period
    duration > half_period && error("duration cannot be longer than half the period")
    b_norm = b * r_star
    speed = 2 * sqrt(r_star^2 - b_norm^2) / duration
    ref_time =  t0 - half_period
    SimpleOrbit(period, t0, b, duration, r_star, b_norm, speed, half_period, ref_time)
end


period(orbit::SimpleOrbit) = orbit.period
duration(orbit::SimpleOrbit) = orbit.duration

planet_position(orbit::SimpleOrbit, t) = relative_position(orbit, t)

function relative_position(orbit::SimpleOrbit, t)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    x = orbit.speed * dt
    y = orbit.b_norm
    z = abs(dt) < 0.5 * orbit.duration ? one(x) : -one(x)
    return SA[x, y, z]
end

function in_transit(orbit::SimpleOrbit, t, r; texp=0)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    tol = sqrt((r + orbit.r_star)^2 - orbit.b_norm^2) / orbit.speed
    tol += 0.5 * texp
    return abs(dt) < tol
end

function in_transit(orbit::SimpleOrbit, t; texp=0)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    tol = 0.5 * (orbit.duration + texp)
    return abs(dt) < tol
end

function flip(orbit::SimpleOrbit, r_planet)
    t0 = orbit.t0 + orbit.half_period
    b = orbit.b_norm / r_planet
    duration = 2 * sqrt(r_planet^2 + orbit.b_norm^2) / orbit.speed
    return SimpleOrbit(orbit.period,
                       t0,
                       b,
                       duration,
                       r_planet,
                       orbit.b_norm,
                       orbit.speed,
                       orbit.half_period,
                       orbit.t0)
end