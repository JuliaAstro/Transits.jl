### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 888afc3a-7715-11eb-3dbe-096fba0a014d
using UnitfulRecipes

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
end

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless example"

# ╔═╡ 22a91084-78a8-11eb-1602-83e6d7cee1bf
KeplerianOrbit(
	ρₛ = 2.37,
	Rₛ = 0.87*6.957e10,
	P =  4.234520*86_400,
	ecc = 0.12,
	t₀ = 0.0,
	incl = 88.60 * π / 180.0,
)

# ╔═╡ 4ba4e3a0-78a3-11eb-13f0-93334f4796a6
KeplerianOrbit(
	ρₛ = 2.37u"g/cm^3",
	Rₛ = 0.87u"Rsun",
	P = 4.234520u"d",
	ecc = 0.12,
	t₀ = 0.0u"d",
	incl = 88.60u"°",
)

# ╔═╡ 01c1cd04-78ab-11eb-2bc9-45e070b29c3d
KeplerianOrbit(
	aRₛ = 13.09,
	P =  4.234520*86_400,
	b = 0.32,
	t₀ = 0.0,
	ecc = 0.12,
)

# ╔═╡ 8a142fca-78aa-11eb-2c6b-112f743d3563
KeplerianOrbit(
	aRₛ = 13.09,
	P = 4.234520u"d",
	b = 0.32, # From separate ρₛ calc. above
	t₀ = 0.0u"d",
	ecc = 0.12,
)

# ╔═╡ accde7a8-76c8-11eb-1584-51b8cf5343e1
let
	orbit = KeplerianOrbit(
		ρₛ = 2.37,
		Rₛ = 0.87*6.957e10,
		P =  4.234520*86_400,
		ecc = 0.12,
		t₀ = 0.0,
		incl = 88.60 * π / 180.0,
	)
	u = [0.4, 0.26] # Quad limb dark
	ld = PolynomialLimbDark(u)
	t = range(-1*86_400, 5*86_400, length=1000) # Days from t0
	rs = range(0, 0.2, length=5) # Radius ratio
	fluxes = @. ld(orbit, t, rs')
	plot(fluxes, label=rs', legend=:bottom, legendtitle=:rprs)
end

# ╔═╡ cf2f7320-752e-11eb-0991-1b852dfc264b
md"""
## "Unit" test: ρₛ vs. aRₛ

Data from: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/hat-p-26)
"""

# ╔═╡ bb362e5a-7719-11eb-10f1-9b877cd2d255
md"### ρₛ"

# ╔═╡ 55b5d50e-76c6-11eb-1e29-39335afb3e59
# Stassun et al. (2017)
HATP26_ρₛ = KeplerianOrbit(
	ρₛ = 2.37u"g/cm^3",
	Rₛ = 0.87u"Rsun",
	P = 4.234520u"d",
	ecc = 0.12,
	t₀ = 0.0u"d",
	incl = 88.60u"°",
)

# ╔═╡ c69d99b6-7714-11eb-1461-4913993fa706
let
	u = [0.4, 0.26]
	ld = PolynomialLimbDark(u)
	t = range(-1, 5, length=1000)u"d"
	rs = range(0, 0.2, length=5)
	fluxes = @. ld(HATP26_ρₛ, t, rs')
	plot(t, fluxes, label=rs', legend=:bottom, legendtitle=:rprs)
end

# ╔═╡ c55a4888-7719-11eb-2e15-bf92d0f117d9
md"### aRₛ"

# ╔═╡ 04de6876-74ff-11eb-18f5-8312aa5218c7
# Stassun et al. (2017)
HATP26_aRs = KeplerianOrbit(
	P = 4.234520u"d",
	aRₛ = 13.09,
	b = 0.32, # From separate ρₛ calc. above
	t₀ = 0.0u"d",
	ecc = 0.12,
)

# ╔═╡ a4d16b3e-7715-11eb-0cb1-c3e2348a87d5
let
	u = [0.4, 0.26]
	ld = PolynomialLimbDark(u)
	t = range(-1, 5, length=1000)u"d"
	rs = range(0, 0.2, length=5)
	fluxes = @. ld(HATP26_aRs, t, rs')
	plot(t, fluxes, label=rs', legend=:bottom, legendtitle=:rprs)
end

# ╔═╡ 9a5edfd2-7795-11eb-39b1-f3a104ab452e
md"## Reduce to simple orbit check"

# ╔═╡ bd9d1a26-777d-11eb-37a3-0defc59efb4a
check_aRₛ(P, δ, T, τ) = (P*δ^(1/4) / (2*π)) * (4 / (T*τ))^(1/2)

# ╔═╡ 72c01acc-7779-11eb-09b3-e7f9e1706b3a
let
	orbit = KeplerianOrbit(
		P = 3.0u"d",
		aRₛ = check_aRₛ(3u"d", 0.2^2, 1u"d", 1u"hr") |> upreferred,
		b = 0.0,
		t₀ = 0.0u"d",
		ecc = 0.0,
	)
	
	u = [0.4, 0.26]
	ld = PolynomialLimbDark(u)
	t = range(-8, 8, length=1000)u"d"
	fluxes = @. ld(orbit, t, 0.2)
	plot(t, fluxes, label=false)
end

# ╔═╡ 0c34711e-7779-11eb-1286-5f9381c67d31
let
	orbit = SimpleOrbit(period=3, duration=1)
	u = [0.4, 0.26]
	ld = PolynomialLimbDark(u)
	t = range(-8, 8, length=1000)
	fluxes = @. ld(orbit, t, 0.2)
	plot(t, fluxes, legend=false, legendtitle=:rprs)
end

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═22a91084-78a8-11eb-1602-83e6d7cee1bf
# ╠═4ba4e3a0-78a3-11eb-13f0-93334f4796a6
# ╠═01c1cd04-78ab-11eb-2bc9-45e070b29c3d
# ╠═8a142fca-78aa-11eb-2c6b-112f743d3563
# ╠═accde7a8-76c8-11eb-1584-51b8cf5343e1
# ╟─cf2f7320-752e-11eb-0991-1b852dfc264b
# ╠═888afc3a-7715-11eb-3dbe-096fba0a014d
# ╟─bb362e5a-7719-11eb-10f1-9b877cd2d255
# ╠═55b5d50e-76c6-11eb-1e29-39335afb3e59
# ╠═c69d99b6-7714-11eb-1461-4913993fa706
# ╟─c55a4888-7719-11eb-2e15-bf92d0f117d9
# ╠═04de6876-74ff-11eb-18f5-8312aa5218c7
# ╠═a4d16b3e-7715-11eb-0cb1-c3e2348a87d5
# ╟─9a5edfd2-7795-11eb-39b1-f3a104ab452e
# ╠═72c01acc-7779-11eb-09b3-e7f9e1706b3a
# ╠═bd9d1a26-777d-11eb-37a3-0defc59efb4a
# ╠═0c34711e-7779-11eb-1286-5f9381c67d31
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
