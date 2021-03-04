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
	Mₛ = 0.151
	Rₛ = 0.189
	P = 0.4626413
	t₀ = 0.2
	b = 0.5
	ecc = 0.1
	ω = 0.1
end

# ╔═╡ 1e7084b4-7cab-11eb-1ab2-97ad504d12fb
orbit = KeplerianOrbit(
	Mₛ = ustrip(u"g", Mₛ * u"Msun"),
	Rₛ = ustrip(u"cm", Rₛ * u"Rsun"),
	P =  ustrip(u"s", P * u"d"),
	t₀ = ustrip(u"s", t₀ * u"d"),
	b = b,
	ecc = ecc,
	ω = ω,
)

# ╔═╡ 415d8962-7d1f-11eb-27c9-5f40fdb3be4c
orbit.incl

# ╔═╡ 2ca5daa6-7d1f-11eb-1715-7d1cd314abce
orbit.ρₛ

# ╔═╡ c8e2fefe-7d1e-11eb-20f7-eb5cb8c21373
uconvert(u"Rsun", orbit.Rₛ * u"cm")

# ╔═╡ e5437dee-7d1e-11eb-2cf3-83fa8723d16c
uconvert(u"d", orbit.t₀ * u"s")

# ╔═╡ f58ef480-7d1e-11eb-023f-13e28437257c
uconvert(u"d", orbit.P * u"s")

# ╔═╡ 1d23c168-7d17-11eb-2cf6-c5ab421b070b
const G_grav = 2942.2062175044193

# ╔═╡ 21dcebe8-7d18-11eb-226f-4f6a20b44f3a
const G = PhysicalConstants.CODATA2018.G

# ╔═╡ 47538ce8-7d17-11eb-0137-8582ad0c6301
a = (G_grav * Mₛ * P^2 / (4 * π^2))^(1.0 / 3) # in Rsun

# ╔═╡ d9ec5c9c-7d17-11eb-2236-5d8674348c95
a_Transits = uconvert(u"cm", (G * 0.151 * u"Msun" * (P * u"d")^2 / (4 * π^2))^(1.0 / 3))

# ╔═╡ 35ff4404-7d24-11eb-2410-03d802c79941
a_Transits

# ╔═╡ af28d4c0-7d19-11eb-2123-316cb6a5a579
uconvert(u"Rsun", a_Transits)

# ╔═╡ 6768fffe-7d1c-11eb-2043-032daafbea9c
uconvert(u"Rsun", orbit.a * u"cm")

# ╔═╡ c494fda4-7d1c-11eb-018d-e13cda2867f5
uconvert(u"Rsun", (orbit.aRₛ * orbit.Rₛ) * u"cm")

# ╔═╡ da3c2c44-7d1d-11eb-1032-117da0ccf171
orbit.aRₛ * 0.189

# ╔═╡ de9fcac6-7d1c-11eb-3f75-13c0465b8864
(orbit.a * u"cm") / (orbit.Rₛ * u"cm")

# ╔═╡ e8cf7b42-7b23-11eb-12b9-9978504654a8
begin
	py"""
	import numpy as np
	from batman import _rsky
	
	def small_star():
		t = np.linspace(0, $P, 500)
	
		a = $(orbit.aRₛ) * $Rₛ # `exoplanet` stores a in solar units
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

# ╔═╡ e5de173a-7d1b-11eb-335d-0b7f19b0fd60
uconvert(u"Rsun", orbit.a * u"cm")

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

# ╔═╡ 722aceea-7cae-11eb-06d7-978d2695c7e9
sum(m) > 0

# ╔═╡ ec86643a-7b23-11eb-0e9f-7df7963f0452
function compute_r(orbit, t)
	pos = relative_position.(orbit, t)
	r = map(pos) do arr
		√(arr[1]^2 + arr[2]^2)
	end
	return r
end

# ╔═╡ 225cdf88-7d22-11eb-048c-a7ae2f92f970
7.0514 * orbit.aRₛ

# ╔═╡ a026e99e-7d1e-11eb-36e8-8d0389fee028
relative_position.(orbit, t)

# ╔═╡ 02dc0966-7cb1-11eb-08d1-e714b6823147
r = compute_r(orbit, t)

# ╔═╡ 13c46c6e-7cb1-11eb-1a0d-413fb58dd2fb
allclose(r_batman[m], r[m], atol=2e-5)

# ╔═╡ 60f3252a-7cb1-11eb-055b-cbbdde15f0a4
r_batman[m]

# ╔═╡ 641aacf0-7cb1-11eb-2fb4-cf81c5d0053c
r[m]

# ╔═╡ 6852e742-7cb1-11eb-291a-13e7eba5b039
plot(r[m])

# ╔═╡ 8ab92060-7cb1-11eb-2a1a-03e986f95f83
plot(r_batman[m])

# ╔═╡ b8c1644c-7d24-11eb-2705-f330a7ad2a32
r[m] ./ r_batman[m] # I will find you >=\

# ╔═╡ 7e4e548c-7d15-11eb-2960-5dfd53827503
plotly()

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═0fc1e9c4-7d18-11eb-3632-ff81656b581d
# ╠═e6012662-7d11-11eb-2986-91b35a94a153
# ╠═5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
# ╠═1e7084b4-7cab-11eb-1ab2-97ad504d12fb
# ╠═415d8962-7d1f-11eb-27c9-5f40fdb3be4c
# ╠═2ca5daa6-7d1f-11eb-1715-7d1cd314abce
# ╠═c8e2fefe-7d1e-11eb-20f7-eb5cb8c21373
# ╠═e5437dee-7d1e-11eb-2cf3-83fa8723d16c
# ╠═f58ef480-7d1e-11eb-023f-13e28437257c
# ╠═1d23c168-7d17-11eb-2cf6-c5ab421b070b
# ╠═21dcebe8-7d18-11eb-226f-4f6a20b44f3a
# ╠═47538ce8-7d17-11eb-0137-8582ad0c6301
# ╠═35ff4404-7d24-11eb-2410-03d802c79941
# ╠═d9ec5c9c-7d17-11eb-2236-5d8674348c95
# ╠═af28d4c0-7d19-11eb-2123-316cb6a5a579
# ╠═6768fffe-7d1c-11eb-2043-032daafbea9c
# ╠═c494fda4-7d1c-11eb-018d-e13cda2867f5
# ╠═da3c2c44-7d1d-11eb-1032-117da0ccf171
# ╠═de9fcac6-7d1c-11eb-3f75-13c0465b8864
# ╠═e8cf7b42-7b23-11eb-12b9-9978504654a8
# ╠═e5de173a-7d1b-11eb-335d-0b7f19b0fd60
# ╠═ef22798c-7cad-11eb-0e44-2b4ea6d7a264
# ╠═9852eaca-7d1f-11eb-013f-c1e8f141d731
# ╠═a216cc32-7d1f-11eb-08f0-efaf94411185
# ╠═4596be2c-7cb1-11eb-364c-ada7577fd31d
# ╠═4441ae7a-7d1e-11eb-3715-2988a36494ed
# ╠═722aceea-7cae-11eb-06d7-978d2695c7e9
# ╠═ec86643a-7b23-11eb-0e9f-7df7963f0452
# ╠═225cdf88-7d22-11eb-048c-a7ae2f92f970
# ╠═a026e99e-7d1e-11eb-36e8-8d0389fee028
# ╠═02dc0966-7cb1-11eb-08d1-e714b6823147
# ╠═13c46c6e-7cb1-11eb-1a0d-413fb58dd2fb
# ╠═60f3252a-7cb1-11eb-055b-cbbdde15f0a4
# ╠═641aacf0-7cb1-11eb-2fb4-cf81c5d0053c
# ╠═6852e742-7cb1-11eb-291a-13e7eba5b039
# ╠═8ab92060-7cb1-11eb-2a1a-03e986f95f83
# ╠═b8c1644c-7d24-11eb-2705-f330a7ad2a32
# ╠═7e4e548c-7d15-11eb-2960-5dfd53827503
