### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	import Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using PlutoUI
	using PyCall
	using Transits.Orbits: KeplerianOrbit, relative_position
	using Unitful, UnitfulAstro
	using StatsPlots
end

# ╔═╡ e34d610a-7b23-11eb-0132-55ce740f8e17
pyimport("pip")["main"](["install", "batman-package", "numpy==1.20.1"])

# ╔═╡ 5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
begin
	Rₛ = 0.189
	P = 0.4626413
	t₀ = 0.2
	b = 0.5
	ecc = 0.1
	ω = 0.1
end

# ╔═╡ 1e7084b4-7cab-11eb-1ab2-97ad504d12fb
orbit = KeplerianOrbit(
	Rₛ = ustrip(u"cm", Rₛ * u"Rsun"),
	P =  ustrip(u"s", P * u"d"),
	t₀ = ustrip(u"s", t₀ * u"d"),
	b = b,
	ecc = ecc,
	ω = ω,
)

# ╔═╡ e8cf7b42-7b23-11eb-12b9-9978504654a8
begin
	py"""
	import numpy as np
	from batman import _rsky
	
	def small_star():
		t = np.linspace(0, $(orbit.P), 500)
	
		a = $(orbit.a)
		incl = $(orbit.incl)

		r_batman = _rsky._rsky(
			t,
			$t₀,
			$P,
			$(orbit.aRₛ),
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

# ╔═╡ 278b5d02-7cb1-11eb-0e3f-d90fe339491e
m = test_vals["m"]

# ╔═╡ 3306b0b4-7cb1-11eb-1d3c-993f49ef2006
t = test_vals["t"]

# ╔═╡ 4596be2c-7cb1-11eb-364c-ada7577fd31d
r_batman = test_vals["r_batman"]

# ╔═╡ 722aceea-7cae-11eb-06d7-978d2695c7e9
m |> sum

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

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═e34d610a-7b23-11eb-0132-55ce740f8e17
# ╠═5fbf61ca-7caf-11eb-19be-7f3cf3e04fa2
# ╠═1e7084b4-7cab-11eb-1ab2-97ad504d12fb
# ╠═e8cf7b42-7b23-11eb-12b9-9978504654a8
# ╠═ef22798c-7cad-11eb-0e44-2b4ea6d7a264
# ╠═278b5d02-7cb1-11eb-0e3f-d90fe339491e
# ╠═3306b0b4-7cb1-11eb-1d3c-993f49ef2006
# ╠═4596be2c-7cb1-11eb-364c-ada7577fd31d
# ╠═722aceea-7cae-11eb-06d7-978d2695c7e9
# ╠═ec86643a-7b23-11eb-0e9f-7df7963f0452
# ╠═02dc0966-7cb1-11eb-08d1-e714b6823147
# ╠═13c46c6e-7cb1-11eb-1a0d-413fb58dd2fb
# ╠═60f3252a-7cb1-11eb-055b-cbbdde15f0a4
# ╠═641aacf0-7cb1-11eb-2fb4-cf81c5d0053c
# ╠═6852e742-7cb1-11eb-291a-13e7eba5b039
# ╠═8ab92060-7cb1-11eb-2a1a-03e986f95f83
