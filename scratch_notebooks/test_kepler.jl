### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	using AstroLib: trueanom, kepler_solver
	using PhysicalConstants
	using PlutoUI
	using PyCall
	using Transits.Orbits: KeplerianOrbit, relative_position, 
	_position, compute_true_anomaly, _star_position, _planet_position, flip
	using Unitful, UnitfulAstro
	using StatsPlots
end

# ╔═╡ 132e3b14-7d78-11eb-3591-bb0f1288b88e
using Test

# ╔═╡ e396ab44-7d2b-11eb-2f91-7573ef462559
md"""
#### Test
"""

# ╔═╡ 8ff6ff7e-7f13-11eb-0b54-83b3ef2ff52b
stack(arr_arr) = hcat((reshape(map(p -> p[i], arr_arr), :) for i in 1:3)...)

# ╔═╡ d2cf7474-7fd2-11eb-1407-151bde70042b
allclose(a, b; rtol=1e-5, atol=1e-8) = all(@. abs(a - b) ≤ (atol + rtol*abs(b)))

# ╔═╡ d86940c0-7ea3-11eb-3082-352f62f765d8
function test_flip()
	t = ustrip.(u"s", range(0, 100; length=1_000)u"d")
	
	orbit = KeplerianOrbit(
		Mₛ = ustrip(u"g", 1.3 * u"Msun"),
		Mₚ = ustrip(u"g", 0.1 * u"Msun"),
		Rₛ = ustrip(u"cm", 1.0 * u"Rsun"),
		P = ustrip(u"s", 100.0 * u"d"),
		t₀ = ustrip(u"s", 0.5 * u"d"),
		incl = 0.25 * π,
		ecc = 0.0,
		ω = 0.5,
		Ω = 1.0
	)
	orbit_flipped = flip(orbit, ustrip(u"cm", 0.7 * u"Rsun"))
		

	u_star = stack(_star_position.(orbit, orbit.Rₛ, t))
	u_planet_flipped = stack(_planet_position.(orbit_flipped, orbit.Rₛ, t))
	for i in 1:3
		@test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
	end
	
	u_planet = stack(_planet_position.(orbit, orbit.Rₛ, t))
	u_star_flipped = stack(_star_position.(orbit_flipped, orbit.Rₛ, t))
	for i in 1:3
		@test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
	end
	
end

# ╔═╡ 6dabd3f2-7fd6-11eb-29d5-bd2620db7b78
test_flip()

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╟─e396ab44-7d2b-11eb-2f91-7573ef462559
# ╠═132e3b14-7d78-11eb-3591-bb0f1288b88e
# ╠═8ff6ff7e-7f13-11eb-0b54-83b3ef2ff52b
# ╠═d86940c0-7ea3-11eb-3082-352f62f765d8
# ╠═d2cf7474-7fd2-11eb-1407-151bde70042b
# ╠═6dabd3f2-7fd6-11eb-29d5-bd2620db7b78
