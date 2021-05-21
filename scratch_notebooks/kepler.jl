### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using PlutoUI
	using PyCall
	using Transits
	using CairoMakie
end

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless examples"

# ╔═╡ da4c42b9-a22b-435d-8c1d-c4b93a9cd411
KeplerianOrbit(
	ρₛ = 2.0,
	Rₛ = 0.5,
	P = 2.0,
	ecc = 0.0,
	t₀ = 0.0,
	incl = π / 2.0,
	Ω = 0.0,
	ω = 0.0,
)

# ╔═╡ 7f40987b-ef60-4786-8639-2100628c4595
KeplerianOrbit(
	rho_s = 2.0,
	Rₛ = 0.5,
	period = 2.0,
	ecc = 0.0,
	t_0 = 0.0,
	incl = π / 2.0,
	Omega = 0.0,
	omega = 0.0,
)

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═7f40987b-ef60-4786-8639-2100628c4595
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
