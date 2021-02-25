### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ fe341026-7278-11eb-2490-0f2ffdeae45e
begin
	using Revise
	using Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using StatsPlots
	using Transits
	using Unitful, UnitfulAstro
end

# ╔═╡ 888afc3a-7715-11eb-3dbe-096fba0a014d
using UnitfulRecipes

# ╔═╡ 4f8e1334-7712-11eb-3da0-f56e7af60908
Unitful.preferunits(u"AU")

# ╔═╡ 785eb51c-76c6-11eb-29fa-e30a76b44588
orbit = KeplerianOrbit(
	ρₛ=2.37,
	Rₛ=0.87*6.957e10,
	P=4.234520*86_400,
	ecc=0.12,
	t₀=0.0,
	incl=88.60 * π / 180,
)

# ╔═╡ accde7a8-76c8-11eb-1584-51b8cf5343e1
u = [0.4, 0.26] # quad limb dark

# ╔═╡ afac1ada-76c8-11eb-327e-539a623394a6
ld = PolynomialLimbDark(u)

# ╔═╡ afac80e2-76c8-11eb-2ddf-a7dae89ed0fb
t = range(-86_400, 5*86_400, length=1000) # seconds from t0

# ╔═╡ afadc18c-76c8-11eb-0aa2-773bc5cdccea
rs = range(0, 0.2, length=10) # radius ratio

# ╔═╡ afb2e158-76c8-11eb-13fd-8112b15851e0
fluxes = @. ld(orbit, t, rs')

# ╔═╡ fa3d19b6-76cb-11eb-0d80-e5e8f2760a6a
plot(fluxes, label=rs', legend=:right, legendtitle=:rprs)

# ╔═╡ cf2f7320-752e-11eb-0991-1b852dfc264b
md"""
#### "Unit" test: ρₛ vs. aRₛ

Data from: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/hat-p-26)
"""

# ╔═╡ bb362e5a-7719-11eb-10f1-9b877cd2d255
md"##### ρₛ"

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
	u = [0.4, 0.26] # quad limb dark
	ld = PolynomialLimbDark(u)
	t = range(-1, 5, length=1000)u"d" # seconds from t0
	rs = range(0, 0.2, length=10) # radius ratio
	fluxes = @. ld(HATP26_ρₛ, t, rs')
	plot(t, fluxes, label=rs', legend=:right, legendtitle=:rprs)
end

# ╔═╡ c55a4888-7719-11eb-2e15-bf92d0f117d9
md"##### aRₛ"

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
	u = [0.4, 0.26] # quad limb dark
	ld = PolynomialLimbDark(u)
	t = range(-1, 5, length=1000)u"d" # seconds from t0
	rs = range(0, 0.2, length=10) # radius ratio
	fluxes = @. ld(HATP26_aRs, t, rs')
	plot(t, fluxes, label=rs', legend=:right, legendtitle=:rprs)
end

# ╔═╡ 9a5edfd2-7795-11eb-39b1-f3a104ab452e
md"##### Simplify to simple orbit"

# ╔═╡ bd9d1a26-777d-11eb-37a3-0defc59efb4a
check_aRₛ(P, δ, T, τ) = (P*δ^(1/4) / (2*π)) * (4 / (T*τ))^(1/2)

# ╔═╡ 72c01acc-7779-11eb-09b3-e7f9e1706b3a
let
	orbit = KeplerianOrbit(
		P = 3u"d",
		aRₛ = check_aRₛ(3u"d", 0.2^2, 1u"d", 1u"hr") |> upreferred,
		b = 0.0,
		t₀ = 0u"d",
		ecc = 0.0,
	)
	#@show orbit
	u = [0.4, 0.26] # quad limb dark
	ld = PolynomialLimbDark(u)
	t = range(-8, 8, length=1000)u"d" # seconds from t0
	fluxes = @. ld(orbit, t, 0.2)
	plot(t, fluxes, label=false)
end

# ╔═╡ 7cec1570-7798-11eb-3328-7d2432ee54a7
KeplerianOrbit(
		P = 3u"d",
		aRₛ = check_aRₛ(3u"d", 0.2^2, 1u"d", 1u"hr") |> upreferred,
		b = 0.0,
		t₀ = 0u"d",
		ecc = 0.0,
	)

# ╔═╡ 0c34711e-7779-11eb-1286-5f9381c67d31
let
	orbit = SimpleOrbit(period=3, duration=1)
	u = [0.4, 0.26] # quad limb dark
	ld = PolynomialLimbDark(u)
	t = range(-8, 8, length=1000)
	fluxes = @. ld(orbit, t, 0.2)
	plot(t, fluxes, label=rs', legend=false, legendtitle=:rprs)
end

# ╔═╡ Cell order:
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
# ╠═4f8e1334-7712-11eb-3da0-f56e7af60908
# ╠═785eb51c-76c6-11eb-29fa-e30a76b44588
# ╠═accde7a8-76c8-11eb-1584-51b8cf5343e1
# ╠═afac1ada-76c8-11eb-327e-539a623394a6
# ╠═afac80e2-76c8-11eb-2ddf-a7dae89ed0fb
# ╠═afadc18c-76c8-11eb-0aa2-773bc5cdccea
# ╠═afb2e158-76c8-11eb-13fd-8112b15851e0
# ╠═fa3d19b6-76cb-11eb-0d80-e5e8f2760a6a
# ╟─cf2f7320-752e-11eb-0991-1b852dfc264b
# ╠═888afc3a-7715-11eb-3dbe-096fba0a014d
# ╟─bb362e5a-7719-11eb-10f1-9b877cd2d255
# ╟─55b5d50e-76c6-11eb-1e29-39335afb3e59
# ╠═c69d99b6-7714-11eb-1461-4913993fa706
# ╟─c55a4888-7719-11eb-2e15-bf92d0f117d9
# ╟─04de6876-74ff-11eb-18f5-8312aa5218c7
# ╠═a4d16b3e-7715-11eb-0cb1-c3e2348a87d5
# ╟─9a5edfd2-7795-11eb-39b1-f3a104ab452e
# ╟─72c01acc-7779-11eb-09b3-e7f9e1706b3a
# ╠═7cec1570-7798-11eb-3328-7d2432ee54a7
# ╟─bd9d1a26-777d-11eb-37a3-0defc59efb4a
# ╠═0c34711e-7779-11eb-1286-5f9381c67d31
