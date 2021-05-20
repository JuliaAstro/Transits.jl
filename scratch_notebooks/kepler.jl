### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ fe341026-7278-11eb-2490-0f2ffdeae45e
begin
	using Revise
	import Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using PlutoUI
	using StatsPlots
	using Transits
	using Unitful
	using UnitfulAstro
	using BenchmarkTools
end

# ╔═╡ 83ed04a9-1913-4b40-afaa-01fa8ecdacd2
using PyCall

# ╔═╡ ac6865fa-47c4-4315-8a4a-398c7909c772
using JLD2

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless examples"

# ╔═╡ e54167b0-55f3-426a-b7ac-a93e8565297f
begin
	make_orbit_ρₛ() = KeplerianOrbit((
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	))
	
	make_orbit_ρₛ_KC() = KeplerianOrbit_KC(
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	)
	
	make_orbit_ρₛ_KD() = KeplerianOrbit_KD(
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	)
	
	orbit_ρₛ,
	orbit_ρₛ_KC,
	orbit_ρₛ_KD = make_orbit_ρₛ(), make_orbit_ρₛ_KC(), make_orbit_ρₛ_KD()
end

# ╔═╡ a6af8efb-ae5e-4f5c-af72-a8687d9bac29
make_orbit_aRₛ() = KeplerianOrbit((
	aRₛ = 7.5,
	P = 2.0,
	incl = π / 2.0,
	# b = 0.0,
	t₀ = 0.0,
	ecc = 0.0,
	Ω = 0.0,
	ω = 0.0,
))

# ╔═╡ f0c2b1b2-b851-4f00-8899-6cf871c190d5
make_orbit_aRₛ()

# ╔═╡ 497322b7-5ff9-4699-8d98-677db5a14061
orbit_aRₛ = make_orbit_aRₛ()

# ╔═╡ b7eb5e36-b4a2-43a8-934d-b76adb794f9f
with_terminal() do
	@btime make_orbit_ρₛ()
	@btime make_orbit_ρₛ_KC()
	@btime make_orbit_ρₛ_KD()
end

# ╔═╡ 5274fb9f-999c-440a-aad2-ce5356157b34
md"""
## Light curves
"""

# ╔═╡ 1a94a407-7a09-4c5f-8796-2570d41ddb4d
t = -0.1:0.01:0.1 # days

# ╔═╡ f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
u = [0.4, 0.26] # Quadratic LD coeffs

# ╔═╡ 638405b7-0f1e-4eb1-85d6-d17dc04f8699
ld = PolynomialLimbDark(u)

# ╔═╡ eb9d7c0c-0790-466c-a722-06c9b41b96cd
begin
	fluxes_ρₛ = @. ld(orbit_ρₛ, t, 0.2)
	fluxes_aRₛ = @. ld(orbit_aRₛ, t, 0.2)
	
	p = plot(xlabel="Time (days)", ylabel="Relative flux")
	plot!(p, t, fluxes_ρₛ, lw=5, label="ρₛ")
	#plot!(p, t, fluxes_aRₛ, label="aRₛ")
end

# ╔═╡ be821fee-6b13-4f02-b7e8-04ed7d56ad39
rotations = [27.014000741644374, 36.63731157445627, 43.95397869463099, 48.503353012894834, 49.99901106593417, 48.34678791380008, 43.650705658060005, 36.20642432789478, 26.482627458951885]

# ╔═╡ 96dc6baf-cc69-4797-a4ac-784feb6897bb
manual = [27.014000741644377, 36.63731157445628, 43.95397869463099, 48.503353012894834, 49.99901106593417, 48.34678791380008, 43.650705658060005, 36.20642432789477, 26.48262745895188]

# ╔═╡ b4dc2635-5a2d-4834-9756-d7a26c007568
separations = 3

# ╔═╡ 6a0e0647-eeb7-491e-b54d-ea8026bad8ed
time = 5:10

# ╔═╡ d9e844b7-d4fa-4ad9-a24f-43da86550560
plot(pos_manual[:, 3] - pos_rotations[:, 3], label="Δz")

# ╔═╡ 5f2c2e95-8962-4b52-8de8-70f186475cb4
begin
	py"""
	import numpy as np
	from batman import _rsky
	
	sky_coords = {}
	def sky_coords():
		t = np.linspace(-100, 100, 1_000)
	
		t0, period, a, e, omega, incl = (
			x.flatten()
			for x in np.meshgrid(
				np.linspace(-5.0, 5.0, 2),
				np.exp(np.linspace(np.log(5.0), np.log(50.0), 3)),
				np.linspace(50.0, 100.0, 2),
				np.linspace(0.0, 0.9, 5),
				np.linspace(-np.pi, np.pi, 3),
				np.arccos(np.linspace(0, 1, 5)[:-1]),
			)
		)
	
		r_batman = np.empty((len(t), len(t0)))
	
		for i in range(len(t0)):
			r_batman[:, i] = _rsky._rsky(
				t, t0[i], period[i], a[i], incl[i], e[i], omega[i], 1, 1
			)
	
		m = r_batman < 100.0
	
		return {
			"m_sum" : m.sum().item(), # Save native Int format
			"r_batman" : r_batman,
			"m" : m,
			"t" : t,
			"t0" : t0,
			"period" : period,
			"a" : a,
			"e" : e,
			"omega" : omega,
			"incl" : incl,
		}
	"""
	sky_coords = py"sky_coords"()
end

# ╔═╡ 02de37b8-d463-44df-9b69-f894f5efa6ee


# ╔═╡ 5e1ceabb-c387-4ed9-ab68-b067a3f1aa92
load("../python_code/test_data/KeplerianOrbit_sky_coords.jld2")

# ╔═╡ 71b217e6-1b19-4e70-b0cd-1ff0637892db
as_matrix(pos) = reinterpret(reshape, Float64, pos) |> permutedims

# ╔═╡ 1785b074-4b48-469f-88d1-0dcc84703f49
pos = Orbits.relative_position.(orbit_ρₛ, t) |> as_matrix

# ╔═╡ 92e87956-d801-426a-a375-3eec4d6f70d2
a, b, c = eachcol(pos)

# ╔═╡ 5b0c03d7-150d-477c-9eae-3df9e46fcd3a
c

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═e54167b0-55f3-426a-b7ac-a93e8565297f
# ╠═a6af8efb-ae5e-4f5c-af72-a8687d9bac29
# ╠═f0c2b1b2-b851-4f00-8899-6cf871c190d5
# ╠═497322b7-5ff9-4699-8d98-677db5a14061
# ╠═b7eb5e36-b4a2-43a8-934d-b76adb794f9f
# ╟─5274fb9f-999c-440a-aad2-ce5356157b34
# ╠═1a94a407-7a09-4c5f-8796-2570d41ddb4d
# ╠═f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
# ╠═638405b7-0f1e-4eb1-85d6-d17dc04f8699
# ╠═eb9d7c0c-0790-466c-a722-06c9b41b96cd
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
# ╠═be821fee-6b13-4f02-b7e8-04ed7d56ad39
# ╠═96dc6baf-cc69-4797-a4ac-784feb6897bb
# ╠═b4dc2635-5a2d-4834-9756-d7a26c007568
# ╠═6a0e0647-eeb7-491e-b54d-ea8026bad8ed
# ╠═d9e844b7-d4fa-4ad9-a24f-43da86550560
# ╠═83ed04a9-1913-4b40-afaa-01fa8ecdacd2
# ╠═5f2c2e95-8962-4b52-8de8-70f186475cb4
# ╠═ac6865fa-47c4-4315-8a4a-398c7909c772
# ╠═02de37b8-d463-44df-9b69-f894f5efa6ee
# ╠═5e1ceabb-c387-4ed9-ab68-b067a3f1aa92
# ╠═1785b074-4b48-469f-88d1-0dcc84703f49
# ╠═92e87956-d801-426a-a375-3eec4d6f70d2
# ╠═5b0c03d7-150d-477c-9eae-3df9e46fcd3a
# ╠═71b217e6-1b19-4e70-b0cd-1ff0637892db
