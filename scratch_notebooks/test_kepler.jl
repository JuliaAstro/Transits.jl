### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	#import Pkg
	#Pkg.activate("..")
	#Pkg.activate("../test")
	using PlutoUI
	using PyCall
	#using Conda
	using Transits.Orbits: KeplerianOrbit, relative_position
	using Unitful, UnitfulAstro
	using StatsPlots
end

# ╔═╡ 0fc1e9c4-7d18-11eb-3632-ff81656b581d
using PhysicalConstants

# ╔═╡ e6012662-7d11-11eb-2986-91b35a94a153
#Conda.add("batman-package"; channel="conda-forge")

# ╔═╡ 5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
begin
	# `exoplanet` input params
	m_star = 0.151
	r_star = 0.189
	period = 0.4626413
	t0 = 0.2
	b = 0.5
	ecc = 0.8
	ω = 0.1
end

# ╔═╡ 1e7084b4-7cab-11eb-1ab2-97ad504d12fb
orbit = KeplerianOrbit(
	Rₛ = ustrip(u"cm", r_star * u"Rsun"),
	Mₛ = ustrip(u"g", m_star * u"Msun"),
	P =  ustrip(u"s", period * u"d"),
	t₀ = ustrip(u"s", t0 * u"d"),
	b = b,
	ecc = ecc,
	ω = ω,
)

# ╔═╡ 02dc0966-7cb1-11eb-08d1-e714b6823147
begin
	py"""
	import numpy as np
	"""
	allclose = py"np.allclose"
end

# ╔═╡ c5c83aba-7d53-11eb-3373-bd6eb798f28a
pos = relative_position.(orbit, orbit.t₀)

# ╔═╡ e396ab44-7d2b-11eb-2f91-7573ef462559
md"""
#### Test
"""

# ╔═╡ dc9afb38-7d53-11eb-217d-67418984f592
allclose((√(pos[1]^2 + pos[2]^2)), orbit.b)

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═0fc1e9c4-7d18-11eb-3632-ff81656b581d
# ╠═e6012662-7d11-11eb-2986-91b35a94a153
# ╠═5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
# ╠═1e7084b4-7cab-11eb-1ab2-97ad504d12fb
# ╠═02dc0966-7cb1-11eb-08d1-e714b6823147
# ╠═c5c83aba-7d53-11eb-3373-bd6eb798f28a
# ╟─e396ab44-7d2b-11eb-2f91-7573ef462559
# ╠═dc9afb38-7d53-11eb-217d-67418984f592
