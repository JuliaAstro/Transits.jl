### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using Transits
	using PlutoUI
	using Transits
	using CairoMakie
	using Measurements
	using Unitful
	using UnitfulAstro
end

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

# ╔═╡ cdda90a8-6f94-4c84-9557-27db687d80e9
KeplerianOrbit(
	ρₛ = 2.0,
	Rₛ = 0.5,
	period = 2.0,
	ecc = 0.0,
	t₀ = 0.0,
	
	incl = 90.,
	Ω = 0.,
	ω = 0.,
)

# ╔═╡ ed6dc004-42a8-472a-a3ea-c7d1abbda868
struct yee{A, B}
	a::A
	b::B
	c::B
	d::B
end

# ╔═╡ fe384a90-d142-4414-a2e5-23ce426af513
yee(1, 2, 3, 4)

# ╔═╡ 7f40987b-ef60-4786-8639-2100628c4595
KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0,
	incl = π / 2.0,
	t₀ = 0.0,
	ecc = 0.0,
	Ω = 0.0,
	ω = 0.0,
)

# ╔═╡ Cell order:
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═cdda90a8-6f94-4c84-9557-27db687d80e9
# ╠═ed6dc004-42a8-472a-a3ea-c7d1abbda868
# ╠═fe384a90-d142-4414-a2e5-23ce426af513
# ╠═7f40987b-ef60-4786-8639-2100628c4595
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
