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
	using Transits
end

# ╔═╡ ff058f64-76cb-11eb-362c-41034d441d4a
using StatsPlots

# ╔═╡ 785eb51c-76c6-11eb-29fa-e30a76b44588
orbit = KeplerianOrbit(
	ρₛ=2.37,
	r_star=0.87*6.957e10,
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
t = range(-1, 1, length=1000) # days from t0

# ╔═╡ afadc18c-76c8-11eb-0aa2-773bc5cdccea
rs = range(0, 0.2, length=10) # radius ratio

# ╔═╡ afb2e158-76c8-11eb-13fd-8112b15851e0
fluxes = @. ld(orbit, t, rs')

# ╔═╡ fa3d19b6-76cb-11eb-0d80-e5e8f2760a6a
plot(fluxes)

# ╔═╡ cf2f7320-752e-11eb-0991-1b852dfc264b
md"""
#### ρₛ vs. aRₛ test #TODO: Units

Data from: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/hat-p-26)
"""

# ╔═╡ 6fac7044-76c6-11eb-172c-0973f14b0fe2
# using Unitful, UnitfulAstro

# ╔═╡ 55b5d50e-76c6-11eb-1e29-39335afb3e59
# HATP26_ρₛ = KeplerianOrbit(
# 	ρₛ=2.37u"g/cm^3",
# 	Rₛ=0.87u"Rsun",
# 	P=4.234520u"d",
# 	ecc=0.12,
# 	t₀=0.0,
# 	incl=88.60u"°",
# )

# ╔═╡ 04de6876-74ff-11eb-18f5-8312aa5218c7
# # Stassun et al. (2017)
# HATP26_aRs = KeplerianOrbit(
# 	P=4.234520u"d",
# 	aRₛ=13.09,
# 	b=0.32, # From separate ρₛ calc. above
# 	t₀=0.0,
# 	ecc=0.12,
# )

# ╔═╡ Cell order:
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
# ╠═ff058f64-76cb-11eb-362c-41034d441d4a
# ╠═785eb51c-76c6-11eb-29fa-e30a76b44588
# ╠═accde7a8-76c8-11eb-1584-51b8cf5343e1
# ╠═afac1ada-76c8-11eb-327e-539a623394a6
# ╠═afac80e2-76c8-11eb-2ddf-a7dae89ed0fb
# ╠═afadc18c-76c8-11eb-0aa2-773bc5cdccea
# ╠═afb2e158-76c8-11eb-13fd-8112b15851e0
# ╠═fa3d19b6-76cb-11eb-0d80-e5e8f2760a6a
# ╟─cf2f7320-752e-11eb-0991-1b852dfc264b
# ╠═6fac7044-76c6-11eb-172c-0973f14b0fe2
# ╠═55b5d50e-76c6-11eb-1e29-39335afb3e59
# ╠═04de6876-74ff-11eb-18f5-8312aa5218c7
