### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ fe341026-7278-11eb-2490-0f2ffdeae45e
begin
	#using Revise
	#using Pkg
	#Pkg.activate(".")
	#Pkg.instantiate()
	using PlutoUI
	using KeywordDispatch
	using Parameters
	#using Transits
end

# ╔═╡ 04491430-7520-11eb-110c-e1b4eee2f974
using Unitful, UnitfulAstro, PhysicalConstants

# ╔═╡ d2c4bc5c-74f5-11eb-3277-9d8d3143f955
#const G = 6.6743e-8 # CGS, CODATA 2018
const G = PhysicalConstants.CODATA2018.G

# ╔═╡ d50493e2-9d8c-4d44-9df7-a9dc71ea18d3
begin	
	get_ρₛ(aRₛ, P) = (3 * π / (G * P^2)) * aRₛ^3
	get_ρₛ(a, P, Rₛ) = get_ρₛ(aRₛ(a, Rₛ), P)
	
	@kwdispatch get_aRₛ()
	@kwmethod get_aRₛ(;ρₛ, P) = cbrt(G * P^2 * ρₛ / (3 * π))
	@kwmethod get_aRₛ(;a, P, Rₛ) = aRₛ(get_ρₛ(a, P, Rₛ), P)
	@kwmethod get_aRₛ(;a, Rₛ) = a / Rₛ
	
	get_a(ρₛ, P, Rₛ) = get_a(get_aRₛ(ρₛ=ρₛ, P=P), Rₛ)
	get_a(aRₛ, Rₛ) = aRₛ * Rₛ
	
	get_b(ρₛ, P, sincos_incl) = get_b(get_aRₛ(ρₛ=ρₛ, P=P), sincos_incl)
	get_b(aRₛ, sincos_incl) = aRₛ * sincos_incl[2]
		
	function get_incl(aRₛ, b, ecc, sincos_ω)
		return acos((b/aRₛ) * (1 + ecc*sincos_ω[1])/(1 - ecc^2))
	end
end

# ╔═╡ 46074bca-5f57-426a-8324-06c5450b711f
begin
	struct KeplerianOrbit
		a
		aRₛ
		b
		ecc
		P
		ρₛ
		Rₛ
		n
		t₀
		incl
		sincos_incl
		Ω
		sincos_Ω
		ω
		sincos_ω
	end
	
	@kwdispatch KeplerianOrbit(;
		Omega => Ω,
		omega => ω,
		aRs => aRₛ,
		rho_s => ρₛ,
		aRs => aRₛ,
		Rs => Rₛ,
		t0 => t₀,
		sincos_Omega => sincos_Ω,
		sincos_omega => sincos_ω,
	)
	
	@kwmethod function KeplerianOrbit(;ρₛ, Rₛ, ecc, P, t₀, incl)
		Ω = π / 2
		ω = 0.0
		sincos_incl = sincos(incl)
		a = get_a(ρₛ, P, Rₛ)
		b = get_b(ρₛ, P, sincos_incl)

		
		return KeplerianOrbit(
			uconvert(u"AU", a),
			get_aRₛ(ρₛ=ρₛ, P=P) |> upreferred,
			upreferred(b) |> upreferred,
			ecc,
			P,
			ρₛ,
			Rₛ,
			2 * π / P,
			t₀,
		    uconvert(u"°", incl),
			sincos(incl),
			Ω,
			sincos(Ω),
			ω,
			sincos(ω)
		)
	end
	
	@kwmethod function KeplerianOrbit(;aRₛ, b, ecc, P, t₀)
		Ω = π / 2
		ω = 0.0
		sincos_ω = sincos(ω)
		incl = get_incl(aRₛ, b, ecc, sincos_ω)
		
		return KeplerianOrbit(
			nothing,
			aRₛ |> upreferred,
			b |> upreferred,
			ecc,
			P,
			nothing,
			nothing,
			2 * π / P,
			t₀,
			uconvert(u"°", incl),
			sincos(incl),
			Ω,
			sincos(Ω),
			ω,
			sincos_ω
		)
	end
end

# ╔═╡ cf2f7320-752e-11eb-0991-1b852dfc264b
md"""
#### Quick test

Data from: [Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/overview/hat-p-26)
"""

# ╔═╡ 75f30152-7504-11eb-0c69-236ff845c2b9
# Stassun et al. (2017)
HATP26_ρₛ = KeplerianOrbit(
	ρₛ=2.37u"g/cm^3",
	Rₛ=0.87u"Rsun",
	P=4.234520u"d",
	ecc=0.12,
	t₀=0.0,
	incl=88.60u"°",
)

# ╔═╡ 04de6876-74ff-11eb-18f5-8312aa5218c7
# Stassun et al. (2017)
HATP26_aRs = KeplerianOrbit(
	P=4.234520u"d",
	aRₛ=13.09,
	b=0.32, # From separate ρₛ calc. above
	t₀=0.0,
	ecc=0.12,
)

# ╔═╡ Cell order:
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
# ╠═d2c4bc5c-74f5-11eb-3277-9d8d3143f955
# ╠═d50493e2-9d8c-4d44-9df7-a9dc71ea18d3
# ╠═46074bca-5f57-426a-8324-06c5450b711f
# ╟─cf2f7320-752e-11eb-0991-1b852dfc264b
# ╠═75f30152-7504-11eb-0c69-236ff845c2b9
# ╠═04de6876-74ff-11eb-18f5-8312aa5218c7
# ╠═04491430-7520-11eb-110c-e1b4eee2f974
