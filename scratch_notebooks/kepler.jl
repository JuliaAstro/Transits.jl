### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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

# ╔═╡ 888afc3a-7715-11eb-3dbe-096fba0a014d
using UnitfulRecipes

# ╔═╡ 13772769-09cc-425b-a5bb-ec3b324319aa
using BenchmarkTools

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless example"

# ╔═╡ 22a91084-78a8-11eb-1602-83e6d7cee1bf
KeplerianOrbit(
	2.37, # ρₛ
	0.87, # Rₛ
	4.234520, # P
	0.12, # ecc
	0.0, # t₀
	88.60 * π / 180.0, # incl
)

# ╔═╡ cf2f7320-752e-11eb-0991-1b852dfc264b
md"""
## "Unit" test: ρₛ vs. aRₛ

Data from: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/hat-p-26)
"""

# ╔═╡ bb362e5a-7719-11eb-10f1-9b877cd2d255
md"### ρₛ"

# ╔═╡ 716e5949-2171-486a-abfc-aafb5cec9f0d
# Stassun et al. (2017)
HATP26_ρₛ = KeplerianOrbit(
	2.37u"g/cm^3", # ρₛ
	0.87u"Rsun", # Rₛ
	4.234520u"d", # P
	0.12, # ecc
	0.0u"d", # t₀
	88.60u"°", # incl
)

# ╔═╡ c55a4888-7719-11eb-2e15-bf92d0f117d9
md"### aRₛ"

# ╔═╡ 04de6876-74ff-11eb-18f5-8312aa5218c7
# Stassun et al. (2017)
HATP26_aRs = KeplerianOrbit(
	13.09, # aRₛ
	4.234520u"d", # P
	0.32, # b
	0.0u"d", # t₀
	0.12, # ecc
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
		check_aRₛ(3u"d", 0.2^2, 1u"d", 1u"hr") |> upreferred, # aRₛ
		3.0u"d", # P
		0.0, # b
		0.0u"d", # t₀
		0.0, # ecc
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

# ╔═╡ 74b35657-0897-418a-a132-d2ecd7917b3c
# make_orbit() = KeplerianOrbit(
# 	2.37, # ρₛ
# 	0.87*6.957e10, # Rₛ
# 	4.234520*86_400, # P
# 	0.12, # ecc
# 	0.0, # t₀
# 	88.60 * π / 180.0, # incl
# )

# ╔═╡ 85292b90-54c8-47b0-94e6-998e2961522a
make_orbit() = KeplerianOrbit(
	2.37, # ρₛ
	0.87*6.957e10, # Rₛ
	4.234520*86_400, # P
	0.12, # ecc
	0.0, # t₀
	88.60 * π / 180.0, # incl
)

# ╔═╡ 038fd025-7354-45c6-bae8-b99c1d1a05ff
make_orbit_units() = KeplerianOrbit(
	2.37u"g/cm^3", # ρₛ
	0.87u"Rsun", # Rₛ
	4.234520u"d", # P
	0.12, # ecc
	0.0u"d", # t₀
	88.60u"°", # incl
)

# ╔═╡ b7eb5e36-b4a2-43a8-934d-b76adb794f9f
with_terminal() do
	@btime make_orbit()
	@btime make_orbit_units()
end

# ╔═╡ ab44e378-8e6d-4e2a-9908-bd409531129a
md"""
## Phase estimates
"""

# ╔═╡ caa432a8-be41-43f5-80df-e149769c9de6
md"""
| Param      | Value | Definition
| ------------------ |:----------- |:-------
| ``N_\text{g}``      |  $(@bind N_g Slider(1:5, show_value=true)) | Number of groups
| ``N_\text{g}^{(\text{exp})}``      | $(@bind N_exp_g Slider(2:20, show_value=true)) | Number of exposures per group
| ``t_\text{start}``      |  $(@bind t_begin Slider(-3:0.5:-1, show_value=true)) hr | Time before mid-transit for first group
| ``t^{(\text{exp})}``   | $(@bind t_exp Slider(1:10:300, show_value=true)) s | Single exposure time
"""

# ╔═╡ 94091f9b-855f-4011-bd52-bb08e1e4f47f
md"""
An initial buffer of $(@bind t_buff Scrubbable(0:30, suffix=" minutes"))
"""

# ╔═╡ dbc088ba-f64e-489d-8101-2e74a88f7658
orbit_HP23 = KeplerianOrbit(
	0.99,              # ρₛ
	1.152,             # Rₛ
	2*1.2128867,       # P
	0.0,               # ecc
	0.0,               # t₀
	85.1 * π / 180.0,  # incl
);

# ╔═╡ bfde28c8-9035-4643-a4bf-a7e059dfef6d
function compute_phases(t, P, t0)
	phase = ((t - t0) / P) % 1
	if phase ≥ 0.5
		phase = phase - 1.0
	end
	return phase
end

# ╔═╡ 275c663f-80dc-4c1b-a443-3784b8023eeb
# https://discourse.julialang.org/t/findnearest-function/4143
function searchsortednearest(a, x)
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

# ╔═╡ bfafb655-14eb-41bf-a619-bb763e773110
begin
	p_HP23 = 0.1113
	u_HP23 = [0.4, 0.26]
	ld_HP23 = PolynomialLimbDark(u_HP23)
	t_HP23 = range(0 - 3 / 24, 0 + 3 / 24, length=1000) # days
	ϕ = compute_phases.(t_HP23, orbit_HP23.P/2.0, orbit_HP23.t₀)
	fluxes_HP23 = @. ld_HP23(orbit_HP23, t_HP23, p_HP23)
	
	########
	# Groups
	########
	# In days
	P_HST = 0.0666667 # 96 min
	t_exp_s = t_exp / 3600 / 24 # 237 s
	
	########
	#N_g = 3
	#N_exp_g = 16
	t_group = N_exp_g * t_exp_s
	
	# Phase points sampled
	t = t_HP23
	t_start = t_begin / 24
	t_end = t_start + N_exp_g * t_exp_s
	t_group_1 = range(t_start, step=t_exp_s, length=N_exp_g)
	
	t_groups = [t_group_1]
	for i in 1:(N_g - 1)
		t_s = t_groups[i][end] + P_HST
	    t_e = t_s + t_group - t_exp_s
	    push!(t_groups, range(t_s, step=t_exp_s, length=N_exp_g))
	end
	
	phase_groups = []
	for t_group in t_groups
	    push!(phase_groups, compute_phases.(t_group, orbit_HP23.P/2.0, orbit_HP23.t₀))
	end
	
	phase_groups = hcat(phase_groups...)
	
	p = plot(xlabel="Planetary phase + 1", ylabel="Normalized flux")
	plot!(p, ϕ .+ 1, fluxes_HP23, label=false)
	scatter!(
		p,
		phase_groups .+ 1,
		#ones(size(phase_groups)[1]),
		fluxes_HP23[searchsortednearest.(Ref(t_HP23), hcat(t_groups...))],
		legend=false,
	)
end

# ╔═╡ 98a55ce2-b9d8-4ec3-bd12-b3b1c8275cee
APT_phase_range = compute_phases.([t_groups[begin][begin] - t_buff/60/24, t_groups[begin][begin+1] + t_buff/60/24], orbit_HP23.P/2.0, orbit_HP23.t₀) .+ 1;

# ╔═╡ 67c1b206-c9e9-4a2d-b471-e8ff6ad49c01
md"""
would correspond to an APT phase of **$(APT_phase_range[1]) - $(APT_phase_range[2])** for the start of the first group 🚀
"""

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═22a91084-78a8-11eb-1602-83e6d7cee1bf
# ╟─cf2f7320-752e-11eb-0991-1b852dfc264b
# ╠═888afc3a-7715-11eb-3dbe-096fba0a014d
# ╟─bb362e5a-7719-11eb-10f1-9b877cd2d255
# ╠═716e5949-2171-486a-abfc-aafb5cec9f0d
# ╟─c55a4888-7719-11eb-2e15-bf92d0f117d9
# ╠═04de6876-74ff-11eb-18f5-8312aa5218c7
# ╠═a4d16b3e-7715-11eb-0cb1-c3e2348a87d5
# ╟─9a5edfd2-7795-11eb-39b1-f3a104ab452e
# ╠═72c01acc-7779-11eb-09b3-e7f9e1706b3a
# ╠═bd9d1a26-777d-11eb-37a3-0defc59efb4a
# ╠═0c34711e-7779-11eb-1286-5f9381c67d31
# ╠═74b35657-0897-418a-a132-d2ecd7917b3c
# ╠═85292b90-54c8-47b0-94e6-998e2961522a
# ╠═038fd025-7354-45c6-bae8-b99c1d1a05ff
# ╠═b7eb5e36-b4a2-43a8-934d-b76adb794f9f
# ╠═13772769-09cc-425b-a5bb-ec3b324319aa
# ╟─ab44e378-8e6d-4e2a-9908-bd409531129a
# ╟─caa432a8-be41-43f5-80df-e149769c9de6
# ╟─bfafb655-14eb-41bf-a619-bb763e773110
# ╟─94091f9b-855f-4011-bd52-bb08e1e4f47f
# ╟─67c1b206-c9e9-4a2d-b471-e8ff6ad49c01
# ╟─98a55ce2-b9d8-4ec3-bd12-b3b1c8275cee
# ╟─dbc088ba-f64e-489d-8101-2e74a88f7658
# ╟─fe341026-7278-11eb-2490-0f2ffdeae45e
# ╟─bfde28c8-9035-4643-a4bf-a7e059dfef6d
# ╟─275c663f-80dc-4c1b-a443-3784b8023eeb
