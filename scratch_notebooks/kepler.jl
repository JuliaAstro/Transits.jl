### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using Transits
	using PlutoUI
	using CairoMakie
	using Measurements
	using Unitful
	using UnitfulAstro
end

# ╔═╡ a94d253f-c5b1-4e2c-80dd-c564cebbf507
using BenchmarkTools

# ╔═╡ da4c42b9-a22b-435d-8c1d-c4b93a9cd411
KeplerianOrbit(
	ρₛ = 2.0,
	Rₛ = 0.5,
	period = 2.0,
	ecc = 0.0,
	t₀ = 0.0,
	incl = π / 2.0,
	Ω = 0.0,
	ω = 0.0,
)

# ╔═╡ 72506d79-5550-4e10-9830-ccf224bb60bd
make_orbit() = KeplerianOrbit(
	ρₛ = 2.0u"g/cm^3",
	Rₛ = 0.5u"Rsun",
	period = 2.0u"d",
	ecc = 0.0,
	t₀ = 0.0u"d",
	incl = 90u"°",
	Ω = 0u"°",
	ω = 0u"°",
)

# ╔═╡ c9cae6cc-f791-4616-b182-41101d2d6b2a
make_orbit_aRₛ() = KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0u"d",
	incl = 90.0u"°",
	t₀ = 0.0u"d",
	ecc = 0.0,
	Ω = 0.0u"°",
	ω = 0.0u"°",
)

# ╔═╡ 2089884b-f2b3-40d3-b983-3fb055199311
make_orbit()

# ╔═╡ cdda90a8-6f94-4c84-9557-27db687d80e9
with_terminal() do
	@btime make_orbit()
	@btime make_orbit_aRₛ()
end

# ╔═╡ 75c3f96d-cb7d-4d81-8cac-b62f4f3a31ea
b =  @benchmark make_orbit()

# ╔═╡ ec41d21b-8e6e-421c-8fb2-a3d41e3ecb1d
b.times |> median

# ╔═╡ f0b6deed-d545-4979-9769-7c5b0b011596
b.allocs

# ╔═╡ 7f40987b-ef60-4786-8639-2100628c4595
KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0u"d",
	incl = 90.0u"°",
	t₀ = 0.0u"d",
	ecc = 0.0,
	Ω = 0.0u"°",
	ω = 0.0u"°",
)

# ╔═╡ Cell order:
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═72506d79-5550-4e10-9830-ccf224bb60bd
# ╠═c9cae6cc-f791-4616-b182-41101d2d6b2a
# ╠═2089884b-f2b3-40d3-b983-3fb055199311
# ╠═cdda90a8-6f94-4c84-9557-27db687d80e9
# ╠═75c3f96d-cb7d-4d81-8cac-b62f4f3a31ea
# ╠═ec41d21b-8e6e-421c-8fb2-a3d41e3ecb1d
# ╠═f0b6deed-d545-4979-9769-7c5b0b011596
# ╠═a94d253f-c5b1-4e2c-80dd-c564cebbf507
# ╠═7f40987b-ef60-4786-8639-2100628c4595
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
