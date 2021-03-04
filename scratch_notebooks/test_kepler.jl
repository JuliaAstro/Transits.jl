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
	ecc = 0.1
	ω = 0.1
end

# ╔═╡ 1e7084b4-7cab-11eb-1ab2-97ad504d12fb
orbit = KeplerianOrbit(
	Mₛ = ustrip(u"g", m_star * u"Msun"),
	Rₛ = ustrip(u"cm", r_star * u"Rsun"),
	P =  ustrip(u"s", period * u"d"),
	t₀ = ustrip(u"s", t0 * u"d"),
	b = b,
	ecc = ecc,
	ω = ω,
)

# ╔═╡ f6f4a05c-7d25-11eb-116d-c1d49a27f8df
orbit2 = KeplerianOrbit(
	ρₛ = 31.53136891307936,
	Rₛ = ustrip(u"cm", r_star * u"Rsun"),
	P =  ustrip(u"s", period * u"d"),
	t₀ = ustrip(u"s", t0 * u"d"),
	b = b,
	ecc = ecc,
	ω = ω,	
)

# ╔═╡ da187a3a-7d2a-11eb-3a6b-f19a5033dc9e
md"""
#### Unit shenanigans
"""

# ╔═╡ fd49c1c6-7d2a-11eb-2b23-2dd939210021
md"""
`exoplanet`
"""

# ╔═╡ 1d23c168-7d17-11eb-2cf6-c5ab421b070b
const G_grav = 2942.2062175044193

# ╔═╡ 47538ce8-7d17-11eb-0137-8582ad0c6301
a = (G_grav * m_star * period^2 / (4 * π^2))^(1.0 / 3) # in Rsun

# ╔═╡ 10af4b46-7d2b-11eb-260d-7f901acb18ae
md"""
`Transits.jl`
"""

# ╔═╡ 21dcebe8-7d18-11eb-226f-4f6a20b44f3a
const G = PhysicalConstants.CODATA2018.G

# ╔═╡ d9ec5c9c-7d17-11eb-2236-5d8674348c95
a_Transits = uconvert(
	u"Rsun",
	(G * m_star * u"Msun" * (period * u"d")^2 / (4 * π^2))^(1.0 / 3)
)

# ╔═╡ e8cf7b42-7b23-11eb-12b9-9978504654a8
begin
	py"""
	import numpy as np
	from batman import _rsky
	
	def small_star():
		t = np.linspace(0, $P, 500)
	
		a = $(orbit.aRₛ) #* $Rₛ # `exoplanet` stores a in solar units
		incl = $(orbit.incl)

		r_batman = _rsky._rsky(
			t,
			$t₀,
			$P,
			a,
			$(orbit.incl),
			$ecc,
			$ω,
			1,
			1
		)

		m = r_batman < 100.0
	
		return {
			"t": t,
			"r_batman": r_batman,
			"m": m,
		}
	"""
	small_star = py"small_star"
	allclose = py"np.allclose"
end

# ╔═╡ ef22798c-7cad-11eb-0e44-2b4ea6d7a264
test_vals = small_star()

# ╔═╡ 9852eaca-7d1f-11eb-013f-c1e8f141d731
t_days = test_vals["t"]

# ╔═╡ a216cc32-7d1f-11eb-08f0-efaf94411185
t = ustrip.(u"s", t_days * u"d")

# ╔═╡ 4596be2c-7cb1-11eb-364c-ada7577fd31d
r_batman = test_vals["r_batman"]

# ╔═╡ 4441ae7a-7d1e-11eb-3715-2988a36494ed
m = test_vals["m"]

# ╔═╡ ec86643a-7b23-11eb-0e9f-7df7963f0452
function compute_r(orbit, t)
	pos = relative_position.(orbit, t)
	r = map(pos) do arr
		√(arr[1]^2 + arr[2]^2)
	end
	return r
end

# ╔═╡ 02dc0966-7cb1-11eb-08d1-e714b6823147
r = compute_r(orbit, t)

# ╔═╡ e396ab44-7d2b-11eb-2f91-7573ef462559
md"""
#### Test
"""

# ╔═╡ 722aceea-7cae-11eb-06d7-978d2695c7e9
sum(m) > 0

# ╔═╡ 13c46c6e-7cb1-11eb-1a0d-413fb58dd2fb
allclose(r_batman[m], r[m], atol=2e-5)

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═0fc1e9c4-7d18-11eb-3632-ff81656b581d
# ╠═e6012662-7d11-11eb-2986-91b35a94a153
# ╠═5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
# ╠═1e7084b4-7cab-11eb-1ab2-97ad504d12fb
# ╠═f6f4a05c-7d25-11eb-116d-c1d49a27f8df
# ╟─da187a3a-7d2a-11eb-3a6b-f19a5033dc9e
# ╟─fd49c1c6-7d2a-11eb-2b23-2dd939210021
# ╠═1d23c168-7d17-11eb-2cf6-c5ab421b070b
# ╠═47538ce8-7d17-11eb-0137-8582ad0c6301
# ╟─10af4b46-7d2b-11eb-260d-7f901acb18ae
# ╠═21dcebe8-7d18-11eb-226f-4f6a20b44f3a
# ╠═d9ec5c9c-7d17-11eb-2236-5d8674348c95
# ╠═e8cf7b42-7b23-11eb-12b9-9978504654a8
# ╠═ef22798c-7cad-11eb-0e44-2b4ea6d7a264
# ╠═9852eaca-7d1f-11eb-013f-c1e8f141d731
# ╠═a216cc32-7d1f-11eb-08f0-efaf94411185
# ╠═4596be2c-7cb1-11eb-364c-ada7577fd31d
# ╠═4441ae7a-7d1e-11eb-3715-2988a36494ed
# ╠═ec86643a-7b23-11eb-0e9f-7df7963f0452
# ╠═02dc0966-7cb1-11eb-08d1-e714b6823147
# ╟─e396ab44-7d2b-11eb-2f91-7573ef462559
# ╠═722aceea-7cae-11eb-06d7-978d2695c7e9
# ╠═13c46c6e-7cb1-11eb-1a0d-413fb58dd2fb
