### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	using AstroLib: trueanom
	using PhysicalConstants
	using PlutoUI
	using PyCall
	using Transits.Orbits: KeplerianOrbit, relative_position
	using Unitful, UnitfulAstro
	using StatsPlots
end

# ╔═╡ 132e3b14-7d78-11eb-3591-bb0f1288b88e
using Test

# ╔═╡ 44c8b5ac-7d6c-11eb-17c5-5f22ffe3e385
Es = [0.0, 2 * π, -226.2, -170.4]

# ╔═╡ c36d1a82-7d6d-11eb-3ee4-25f361b42c70
eccs = fill(1.0 - 1e-6, length(Es))

# ╔═╡ c36d8506-7d6d-11eb-3d00-4b0ec4f86707
eccs[end] = 0.9939879759519037

# ╔═╡ 8ed6227a-7d75-11eb-39be-a9d675783c06
eccs

# ╔═╡ 4f8b7e3e-7d7b-11eb-3e36-8b2d52ca6b9a
function compute_M_and_ν(E, ecc)
	M = E - ecc * sin(E)
	ν = trueanom(E, ecc)
	return M, ν
end

# ╔═╡ 3893efb6-7d6e-11eb-32bb-1bfcc375b17d
# E: Eccentric anomaly
# ecc: eccentricity
# Compute sin_ν₀, cos_ν₀ in a numerically stable manner
function compute_sincos_ν_stable(E, ecc; tol=1e-10) # <- E₉
	sin_E, cos_E = sincos(E)

	# Compute cos(ν), sin(ν) where valid
	denom = 1.0 + cos_E
	m = denom > tol
	m_inv = !m
	denom += 1.0 * m_inv

	# First, compute tan(0.5*E) = sin(E) / (1 + cos(E))
	tan_ν₂ = √((1 + ecc) / (1 - ecc)) * sin_E / denom  # tan(0.5*ν)
	tan²_ν₂ = tan_ν₂ * tan_ν₂

	# Then we compute sin(ν) and cos(ν) using:
	#  sin(ν) = 2*tan(0.5*ν)/(1 + tan(0.5*ν)²), and
	#  cos(ν) = (1 - tan(0.5*ν)²)/(1 + tan(0.5*ν)²)
	denom = 1.0 / (1.0 + tan²_ν₂)
	sin_ν = 2.0 * tan_ν₂ * denom
	cos_ν = (1.0 - tan²_ν₂) * denom

	return sin_ν * m, cos_ν * m - 1.0 * m_inv
end

# ╔═╡ e396ab44-7d2b-11eb-2f91-7573ef462559
md"""
#### Test
"""

# ╔═╡ b4ed0446-7d76-11eb-3506-4552b1f0f11a
begin
	py"""
	import numpy as np
	"""
	const allclose = py"np.allclose"
end

# ╔═╡ a7011f58-7d96-11eb-24e0-d5942c10bf88
function kepler_solver(M::T, e::T) where {T<:AbstractFloat}
	@assert 0 <= e <= 1 "eccentricity must be in the range [0, 1]"
    # M must be in the range [-pi, pi], see Markley (1995), page 2.
    #M = rem2pi(M, RoundNearest)
	M = M % (2.0*π) == 0.0 ? 0.0 : mod2pi(M)
    if iszero(M) || iszero(e)
        return M
    else
		high = M > π
		high && (M = 2.0*π - M)
        pi2 = abs2(T(pi))
        α = (3 * pi2 + 1.6 * (pi2 - pi * abs(M))/(1 + e))/(pi2 - 6)
        d = 3 * (1 - e) + α * e
        q = 2 * α * d * (1 - e) - M * M
        r = 3 * α * d * (d - 1 + e) * M + M * M * M
        w = cbrt(abs2(abs(r) + sqrt(q * q * q + r * r)))
        E1 = (2 * r * w / @evalpoly(w, q * q, q, 1) + M)/d
        f2, f3 = e .* sincos(E1)
        f0 = E1 - f2 - M
        f1 = 1 - f3
        δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        δ4 = -f0 / @evalpoly(δ3, f1, f2 / 2, f3 / 6)
        δ5 = -f0 / @evalpoly(δ4, f1, f2 / 2, f3 / 6, - f2 / 24)
        E = E1 + δ5 # equation 29
		high && (E = 2.0*π - E)
		return E
    end
end

# ╔═╡ f6705594-7d6c-11eb-3624-53324bb04681
function compute_E₀(E, ecc)
	#  2 atan(√((1 + e) / (1 - e)) * tan(E / 2)
	M, ν = compute_M_and_ν(E, ecc)
	sin_ν, cos_ν = sincos(ν)
	
	return kepler_solver(M, ecc), sin_ν, cos_ν # Markley (1995)
end

# ╔═╡ 75d06d62-7d73-11eb-3db7-97c3cbaebdb9
function edge(E, ecc)
	E₀, sin_ν, cos_ν = compute_E₀(E, ecc)
	sin_ν₀, cos_ν₀ = compute_sincos_ν_stable(E₀, ecc)
	E_mod = E % (2.0*π) == 0.0 ? 0.0 : mod2pi(E)
	return [E_mod, E₀, sin_ν, sin_ν₀, cos_ν, cos_ν₀]
end

# ╔═╡ 31077a28-7d7e-11eb-340c-77db01116fc1
(
	E_calc,     E_test,
	sin_ν_calc, sin_ν_test,
	cos_ν_calc, cos_ν_test,
) = eachrow(hcat(edge.(Es, eccs)...))

# ╔═╡ a586a76c-7d87-11eb-3d5a-532dabe16cd0
all(isfinite.(sin_ν_test))

# ╔═╡ b365405a-7d87-11eb-3d09-2ff43daad608
all(isfinite.(cos_ν_test))

# ╔═╡ db4c082a-7d7c-11eb-366f-6da594876c0c
E_calc, E_test, allclose(E_calc, E_test)

# ╔═╡ 715065fc-7d80-11eb-0be1-adc0b76315db
sin_ν_calc, sin_ν_test, allclose(sin_ν_calc, sin_ν_test)

# ╔═╡ 762dcf88-7d80-11eb-39eb-fb295bf968d1
cos_ν_calc, cos_ν_test, allclose(cos_ν_calc, cos_ν_test)

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═44c8b5ac-7d6c-11eb-17c5-5f22ffe3e385
# ╠═c36d1a82-7d6d-11eb-3ee4-25f361b42c70
# ╠═c36d8506-7d6d-11eb-3d00-4b0ec4f86707
# ╠═8ed6227a-7d75-11eb-39be-a9d675783c06
# ╠═4f8b7e3e-7d7b-11eb-3e36-8b2d52ca6b9a
# ╠═f6705594-7d6c-11eb-3624-53324bb04681
# ╠═3893efb6-7d6e-11eb-32bb-1bfcc375b17d
# ╟─e396ab44-7d2b-11eb-2f91-7573ef462559
# ╠═132e3b14-7d78-11eb-3591-bb0f1288b88e
# ╠═75d06d62-7d73-11eb-3db7-97c3cbaebdb9
# ╠═31077a28-7d7e-11eb-340c-77db01116fc1
# ╠═a586a76c-7d87-11eb-3d5a-532dabe16cd0
# ╠═b365405a-7d87-11eb-3d09-2ff43daad608
# ╠═db4c082a-7d7c-11eb-366f-6da594876c0c
# ╠═715065fc-7d80-11eb-0be1-adc0b76315db
# ╠═762dcf88-7d80-11eb-39eb-fb295bf968d1
# ╠═b4ed0446-7d76-11eb-3506-4552b1f0f11a
# ╠═a7011f58-7d96-11eb-24e0-d5942c10bf88
