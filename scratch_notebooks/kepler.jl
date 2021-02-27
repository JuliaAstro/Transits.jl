### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 81592330-78a3-11eb-0bb7-29b2ee77a25c
using KeywordDispatch, PhysicalConstants

# ╔═╡ 888afc3a-7715-11eb-3dbe-096fba0a014d
using UnitfulRecipes

# ╔═╡ fe341026-7278-11eb-2490-0f2ffdeae45e
begin
	using Revise
	using Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using PlutoUI
	using StatsPlots
	using Transits
	using Unitful, UnitfulAstro
end

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless example"

# ╔═╡ a4c873ca-78a3-11eb-1110-1f3970469670
begin
	const G = PhysicalConstants.CODATA2018.G
	const G_cgs = ustrip(u"cm^3/g/s^2", G)
	const M₀ = π / 2.0
end

# ╔═╡ b7480a4c-78a3-11eb-1427-5dc77e9b39a3
begin
	compute_ρₛ(aRₛ, P) = (3.0 * π / (G_cgs * P^2.0)) * aRₛ^3.0
	compute_ρₛ(Rₛ, P::T) where {T <: Unitful.Time} = (3.0 * π / (G * P^2.0)) * aRₛ^3.0
	compute_ρₛ(a, P, Rₛ) = compute_ρₛ(aRₛ(a, Rₛ), P)
	
	# Semi-major axis / star radius ratio
	@kwdispatch compute_aRₛ()
	@kwmethod compute_aRₛ(;ρₛ, P) = cbrt(G_cgs * P^2.0 * ρₛ / (3.0 * π))
	@kwmethod compute_aRₛ(;ρₛ, P::T) where {T <: Unitful.Time} = cbrt(G * P^2.0 * ρₛ / (3.0 * π))
	@kwmethod compute_aRₛ(;a, P, Rₛ) = aRₛ(compute_ρₛ(a, P, Rₛ), P)
	@kwmethod compute_aRₛ(;a, Rₛ) = a / Rₛ
	
	# Semi-major axis
	compute_a(ρₛ, P, Rₛ) = compute_a(compute_aRₛ(ρₛ=ρₛ, P=P), Rₛ)
	compute_a(aRₛ, Rₛ) = aRₛ * Rₛ
	
	# Impact parameter
	compute_b(ρₛ, P, sincosi) = compute_b(compute_aRₛ(ρₛ=ρₛ, P=P), sincosi)
	compute_b(aRₛ, sincosi) = aRₛ * sincosi[2]
	
	# Inclination
	function compute_incl(aRₛ, b, ecc, sincosω)
	    return acos((b/aRₛ) * (1.0 + ecc*sincosω[1])/(1.0 - ecc^2))
	end
end

# ╔═╡ 473b2f04-78a3-11eb-3e06-7d44497b7962
begin
	struct Kep{T,L,D,R,A,I}
	    a::L
	    aRₛ::R
	    b::R
	    ecc::R
	    P::T
	    ρₛ::D
	    Rₛ::L
	    n::I
	    t₀::T
	    incl::A
	    Ω::R
	    ω::R
	end
	
	# Enable keyword dispatch and argument name aliasing
	@kwdispatch Kep(;
	    Omega => Ω,
	    omega => ω,
	    aRs => aRₛ,
	    rho_s => ρₛ,
	    aRs => aRₛ,
	    Rs => Rₛ,
	    t0 => t₀,
	)
	
	@kwmethod function Kep(;ρₛ, Rₛ, ecc, P, t₀, incl)
		#Rₛ isa Real && (Rₛ = Rₛ * 6.957e10)
    	#P isa Real && (P = P * 86_400.0)
    	#incl isa Real && (incl = incl * π / 180.0)
	    Ω = M₀
	    ω = 0.0
	    aRₛ = compute_aRₛ(ρₛ=ρₛ, P=P)
	    a = compute_a(ρₛ, P, Rₛ)
	    b = compute_b(ρₛ, P, sincos(incl))
		n = 2.0 * π / P
	
	    # Normalize quantities
		@show ρₛ, P, Rₛ
		a, Rₛ = promote(a, Rₛ)
		
		# Normalize unitless types
	    aRₛ, b, ecc = promote(aRₛ, b, ecc)
		
		#@show typeof(a) typeof(aRₛ) typeof(b) typeof(ecc) typeof(P) typeof(ρₛ) typeof(Rₛ) typeof(n) typeof(t₀) typeof(incl) typeof(Ω) typeof(ω)
			
	    @show a
		return Kep(
	        a,
	        aRₛ,
	        b,
	        ecc,
	        P,
	        ρₛ,
	        Rₛ,
	        n,
	        t₀,
	        incl,
	        Ω,
	        ω,
	    )
	end
	
	@kwmethod function Kep(;aRₛ, P, b, t₀, ecc)
    Ω = M₀
    ω = 0.0
    incl = compute_incl(aRₛ, b, ecc, sincos(ω))

    return Kep(
        nothing,
        aRₛ,
        b,
        ecc,
        P,
        nothing,
        nothing,
        2.0 * π / P,
        t₀,
        incl,
        Ω,
        ω,
    )
	end

end

# ╔═╡ 22a91084-78a8-11eb-1602-83e6d7cee1bf
Kep(
	ρₛ = 2.37,
	Rₛ = 0.87*6.957e10,
	P =  4.234520*86_400,
	ecc = 0.12,
	t₀ = 0.0,
	incl = 88.60 * π / 180.0,
)

# ╔═╡ 4ba4e3a0-78a3-11eb-13f0-93334f4796a6
Kep(
	ρₛ = 2.37u"g/cm^3",
	Rₛ = 0.87u"Rsun",
	P = 4.234520u"d",
	ecc = 0.12,
	t₀ = 0.0u"d",
	incl = 88.60u"°",
)

# ╔═╡ 01c1cd04-78ab-11eb-2bc9-45e070b29c3d
Kep(
	aRₛ = 13.09,
	P =  4.234520*86_400,
	b = 0.32,
	t₀ = 0.0,
	ecc = 0.12,
)

# ╔═╡ 8a142fca-78aa-11eb-2c6b-112f743d3563
Kep(
	aRₛ = 13.09,
	P = 4.234520u"d",
	b = 0.32, # From separate ρₛ calc. above
	t₀ = 0.0u"d",
	ecc = 0.12,
)

# ╔═╡ 785eb51c-76c6-11eb-29fa-e30a76b44588
orbit = KeplerianOrbit(
	ρₛ = 2.37,
	Rₛ = 0.87*6.957e10,
	P =  4.234520*86_400,
	ecc = 0.12,
	t₀ = 0.0,
	incl = 88.60 * π / 180.0,
)

# ╔═╡ accde7a8-76c8-11eb-1584-51b8cf5343e1
let
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
		P = 3u"d",
		aRₛ = check_aRₛ(3u"d", 0.2^2, 1u"d", 1u"hr") |> upreferred,
		b = 0.0,
		t₀ = 0u"d",
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
# ╠═81592330-78a3-11eb-0bb7-29b2ee77a25c
# ╠═a4c873ca-78a3-11eb-1110-1f3970469670
# ╠═473b2f04-78a3-11eb-3e06-7d44497b7962
# ╠═22a91084-78a8-11eb-1602-83e6d7cee1bf
# ╠═4ba4e3a0-78a3-11eb-13f0-93334f4796a6
# ╠═01c1cd04-78ab-11eb-2bc9-45e070b29c3d
# ╠═8a142fca-78aa-11eb-2c6b-112f743d3563
# ╠═b7480a4c-78a3-11eb-1427-5dc77e9b39a3
# ╠═785eb51c-76c6-11eb-29fa-e30a76b44588
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
