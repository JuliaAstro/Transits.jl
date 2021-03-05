### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	using AstroLib: kepler_solver, trueanom
	using PhysicalConstants
	using PlutoUI
	using PyCall
	using Transits.Orbits: KeplerianOrbit, relative_position
	using Unitful, UnitfulAstro
	using StatsPlots
end

# ╔═╡ 44c8b5ac-7d6c-11eb-17c5-5f22ffe3e385
E = [0.0, 2 * π, -226.2, -170.4]

# ╔═╡ c36d1a82-7d6d-11eb-3ee4-25f361b42c70
ecc = fill(1.0 - 1e-6, length(E))

# ╔═╡ c36d8506-7d6d-11eb-3d00-4b0ec4f86707
ecc[end] = 0.9939879759519037

# ╔═╡ f6705594-7d6c-11eb-3624-53324bb04681
function compute_E₀(E, ecc)
	#  2 atan(√((1 + e) / (1 - e)) * tan(E / 2)
	M = E - ecc * sin(E)
	ν = trueanom(E, ecc)
	sin_ν, cos_ν = sincos(ν)
	
	return kepler_solver.(M, ecc), sin_ν, cos_ν # Markley (1995)
end

# ╔═╡ 0a8c32d2-7d6d-11eb-3e1e-adf1d2170ae6
compute_E₀.(E, ecc) # E₀, sin_ν, cos_ν

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
	cos_ν = (1 - tan²_ν₂) * denom
	sin_ν = 2 * tan_ν₂ * denom

	return sin_ν * m, cos_ν * m - 1.0 * m_inv
end

# ╔═╡ d9a077c4-7d72-11eb-3589-8d016fa8d5bb
compute_sincos_ν_stable.(E, ecc)

# ╔═╡ e396ab44-7d2b-11eb-2f91-7573ef462559
md"""
#### Test
"""

# ╔═╡ 75d06d62-7d73-11eb-3db7-97c3cbaebdb9
function edge(E, ecc)
	E₀, sin_ν, cos_ν = compute_E₀(E, ecc) # E₀, sin_ν, cos_ν
	sin_ν₀, cos_ν₀ = compute_sincos_ν_stable(E₀, ecc)
	return [
		E₀         sin_ν₀ cos_ν₀
		mod2pi(E)  sin_ν  cos_ν
	]
end

# ╔═╡ e1230202-7d73-11eb-077b-2bd6aba9a40e
edge.(E, ecc)

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═44c8b5ac-7d6c-11eb-17c5-5f22ffe3e385
# ╠═c36d1a82-7d6d-11eb-3ee4-25f361b42c70
# ╠═c36d8506-7d6d-11eb-3d00-4b0ec4f86707
# ╠═f6705594-7d6c-11eb-3624-53324bb04681
# ╠═0a8c32d2-7d6d-11eb-3e1e-adf1d2170ae6
# ╠═d9a077c4-7d72-11eb-3589-8d016fa8d5bb
# ╠═3893efb6-7d6e-11eb-32bb-1bfcc375b17d
# ╠═e396ab44-7d2b-11eb-2f91-7573ef462559
# ╠═75d06d62-7d73-11eb-3db7-97c3cbaebdb9
# ╠═e1230202-7d73-11eb-077b-2bd6aba9a40e
