

struct SecondaryLimbDark{LD1<:AbstractLimbDark,LD2<:AbstractLimbDark,T} <: AbstractLimbDark
    primary_driver::LD1
    secondary_driver::LD2
    brightness_ratio::T
end

@doc raw"""
"""
function SecondaryLimbDark(driver1::AbstractLimbDark, driver2::AbstractLimbDark; brightness_ratio)
    return SecondaryLimbDark(driver1, driver2, brightness_ratio)
end

function SecondaryLimbDark(u1::AbstractVector, u2=u1; kwargs...)
    driver1 = PolynomialLimbDark(u1)
    driver2 = PolynomialLimbDark(u2)
    return SecondaryLimbDark(driver1, driver2; kwargs...)
end

function compute(ld::SecondaryLimbDark, orbit::AbstractOrbit, t, r)
    f1 = ld.primary_driver(orbit, t, r)
    orbit2 = Orbits.flip(orbit, r * orbit.r_star)
    f2 = ld.secondary_driver(orbit2, t, 1/r)
    flux_ratio = ld.brightness_ratio * r^2
    return (f1 + flux_ratio * f2) / (1 + flux_ratio)
end
