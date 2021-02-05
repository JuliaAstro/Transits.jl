struct SimpleOrbit
    period
    t0
    b
    duration
    Rstar

    b_norm
    speed
    half_period
    ref_time
end

"""
    SimpleOrbit(; period, duration, t0=0, b=0, Rstar=1)

Circular orbit parameterized by the basic observables of a transiting system.

# Parameters
* `period` - The orbital period of the planets, nominally in days
* `duration` The duration of the transit, same units as `period`
* `t0` - The midpoint time of the reference transit, same units as `period`
* `b` - The impact parameter of the orbit
* `Rstar` - The radius of the star, nominally in solar radii
"""
function SimpleOrbit(;period, duration, t0=0.0, b=0.0, Rstar=1.0)
    b_norm = b * Rstar
    speed = 2 * sqrt(Rstar^2 - b_norm^2)
    half_period = period / 2
    ref_time =  t0 - half_period
    SimpleOrbit(period, t0, b, duration, Rstar, b_norm, speed, half_period, ref_time)
end


period(orbit::SimpleOrbit) = orbit.period
duration(orbit::SimpleOrbit) = orbit.duration

planet_position(orbit::SimpleOrbit, t) = relative_position(orbit, t)

function relative_position(orbit::SimpleOrbit, t)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    x = orbit.speed * dt
    y = orbit.b_norm
    z = abs(dt) < 0.5 * duration(orbit) ? one(x) : -one(x)
    return SA[x, y, z]
end

function in_transit(orbit::SimpleOrbit, t, r; texp=0)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    tol = sqrt((r + orbit.Rstar)^2 - orbit.b_norm^2) / orbit.speed
    tol += 0.5 * texp
    return abs(dt) < tol
end

function in_transit(orbit::SimpleOrbit, t; texp=0)
    dt = (t - orbit.ref_time) % period(orbit) - orbit.half_period
    tol = 0.5 * duration(orbit)
    tol += 0.5 * texp
    return abs(dt) < tol
end
