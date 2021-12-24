compute_R_star(rho_star, M_star) = cbrt(3.0 * M_star / (4.0 * π * rho_star))
compute_R_star_nom(G::Real) = 1.0
compute_R_star_nom(G) = 1.0u"Rsun"
compute_M_star(rho_star, R_star) = 4.0 * π * R_star^3 * rho_star / 3.0
compute_M_planet_nom(G::Real) = 0.0
compute_M_planet_nom(G) = 0.0u"Msun"
compute_M_tot(m1, m2) = m1 + m2
compute_M_tot(a, G, period) = 4.0 * π^2 * a^3 / (G * period^2)
compute_a(aR_star, R_star) = R_star * aR_star
compute_a(M_tot, period, G) = cbrt(G * M_tot * period^2 / (4.0 * π^2))
compute_a_X(a, m, M) = a * m / M
compute_period(M_tot, a, G) = 2.0 * π * sqrt(a^3 / (G * M_tot))
compute_rho_star(M_star, R_star) = 3.0 * M_star / (4.0 * π * R_star^3)
function compute_E0(ecc, cos_omega, sin_omega)
    y = sqrt(1.0 - ecc) * cos_omega
    x = sqrt(1.0 + ecc) * (1.0 + sin_omega)
    return 2.0 * atan(y, x)
end
compute_M0(ecc, E0) = E0 - ecc * sin(E0)
compute_M(t, t0, t_ref, n) = (t - t0 - t_ref) * n

# Spherical density
compute_rho(M, R) = 0.75 * M / (π * R^3)
compute_rho(M, R::Nothing) = nothing

# Semi-major axis / star radius ratio, assuming circular orbit
function compute_aor(duration, period, b; r=nothing)
    r = isnothing(r) ? 0.0 : r
    sin_ϕ, cos_ϕ = sincos(π * duration / period)
    return sqrt((1 + r)^2 - (b * cos_ϕ)^2) / sin_ϕ
end

# Impact radius
function compute_b(a_planet, R_star, duration, period, incl_factor_inv, ecc, sin_omega)
    c = sin(π * duration / (period * incl_factor_inv))
    c_sq = c^2
    ecc_sin_omega = ecc*sin_omega
    aor = a_planet / R_star
    num = aor^2 * c_sq - 1.0
    den = c_sq * ecc_sin_omega^2 + 2.0 * c_sq * ecc_sin_omega + c_sq - ecc^4 + 2.0 * ecc^2 - 1.0
    return sqrt(num / den) * (1.0 - ecc) * (1.0 + ecc)
end
compute_b(cos_incl, dcosi_db) = cos_incl / dcosi_db

# Inclination factor
compute_incl_factor_inv(ecc, sin_omega) = (1.0 - ecc) * (1.0 + ecc) / (1.0 + ecc * sin_omega)

# Jacobian for cos(i) -> b
compute_dcosi_db(a, R_star, incl_factor_inv) = R_star / (a * incl_factor_inv)

# Planet radius
compute_R_planet(R_star, r, R_planet) = R_planet
compute_R_planet(R_star, r, R_planet::Nothing) = iszero(r) ? zero(R_star) : R_star * r

# Transit times
compute_t0_tp(t0::Nothing, tp; M0, n) = (tp + M0 / n, tp)
compute_t0_tp(t0, tp::Nothing; M0, n) = (t0, t0 - M0 / n)
compute_t0_tp(t0::Nothing, tp::Nothing; kwargs...) = throw(ArgumentError("Please specify either `t0` or `tp`"))
compute_t0_tp(t0, tp; kwargs...) = throw(ArgumentError("Please only specify one of `t0` or `tp`"))
