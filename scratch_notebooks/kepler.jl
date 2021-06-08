### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using Transits
	using BenchmarkTools
	import PlutoUI as Pl
	using CairoMakie
	using Measurements
	using Unitful
	using UnitfulAstro
end

# ╔═╡ 2a67f029-8472-47f4-8fdc-054214861355
using KeywordCalls

# ╔═╡ afb5a050-b80d-4205-8583-6515779c18e3
KeplerianOrbit(
	ρₛ = 2.0u"g/cm^3",
	Rₛ = 1.0u"Rsun",
	P = 2.0u"d",
	ecc = 0.001,
	t₀ = 0.0u"d",
	incl = 0.0u"°",
	Ω = 0.0u"°",
	ω = 0.0u"°",
)

# ╔═╡ 0895d026-491f-4168-8ad0-509a3227ca89
KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0u"d",
	incl = 90.0u"°",
	t₀ = 0.0u"d",
	ecc = 0.0,	
	Ω = 0.0u"°",
	ω = 0.0u"°",
	#R_s = 2.0u"Rsun",
)

# ╔═╡ da4c42b9-a22b-435d-8c1d-c4b93a9cd411
function make_orbit_ρₛ(P, δ, T, τ)
	return KeplerianOrbit(
	ρₛ = (3.0*P / (Transits.Orbits.G_nom*π^2)) * (δ^(1/4) / (√(T*τ)))^3,
	Rₛ = 1.0,
	P = P,
	ecc = 0.001,
	t₀ = 0.0,
	#b = √(1 - √(δ)*T/τ),
	incl = 0.,
	Ω = 0.0,
	ω = 0.0,
)
end

# ╔═╡ 09d609e7-54c4-468c-8f8f-39aa808d20be
function allclose(a, b; rtol = 1e-5, atol = 1e-8)
    return all(abs.(a - b) .<= (atol .+ rtol * abs.(b)))
end

# ╔═╡ 3de4e3b5-3d64-4e50-b951-de03a8d77ad1
const flip = Orbits.flip

# ╔═╡ 72174e0c-0116-4d1e-8365-c1aa63b648ca
const _star_position = Orbits._star_position

# ╔═╡ 6de491a0-525d-4283-999f-dff478185214
const _planet_position = Orbits._planet_position

# ╔═╡ df13cf06-bb11-452f-978e-98d27e226bb2
as_matrix(pos) = reinterpret(reshape, Float64, pos) |> permutedims

# ╔═╡ 81d987d3-b0e8-48cb-9a9c-d924cd8f7297
let
    t = range(0, 100; length=1_000)
    m_star = 1.3
	r_star = 1.1
    orbit = KeplerianOrbit(
        rho_s = m_star / ((4.0/3.0) * π * r_star^3),
        R_s = r_star,
        period = 100.0,
        ecc = 0.3,
        t_0 = 0.5,
        incl = 0.25*π,
        omega = 0.5,
        Omega = 1.0,
		M_p = 0.1,
    )
    
	orbit_flipped = flip(orbit, 0.7)

    u_star = as_matrix(_star_position.(orbit, orbit.R_s, t))
    u_planet_flipped = as_matrix(_planet_position.(orbit_flipped, orbit.R_s, t))
#     for i in 1:3
#         #@test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
#     end

#     u_planet = stack(_planet_position.(orbit, orbit.R_s, t))
#     u_star_flipped = stack(_star_position.(orbit_flipped, orbit.R_s, t))
#     for i in 1:3
#         #@test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
#     end
	
	#u_star, u_planet_flipped
	
	fig = Figure(resolution=(1_000, 400))
	ax1 = Axis(fig[1, 1])
	ax2 = Axis(fig[1, 2])
	linkyaxes!(ax1, ax2)
	
	for (star1, star2, pos) in zip(eachcol(u_star), eachcol(u_planet_flipped), ["x", "y", "z"])
		lines!(ax1, t, star1, label=pos)
		lines!(ax2, t, star2, label=pos)
	end
	
	axislegend(ax2)
	
	fig
end

# ╔═╡ ca156253-c53d-4b5c-a8b8-7fe694b8f095
h(nt::NamedTuple{(:b, :a, :c)}) = println("Calling f(b = ", nt.b,",a = ", nt.a, ")")

# ╔═╡ 9d2d78b1-2319-43d2-8a51-c383ba17b585
@kwcall h(b,a,c=0)

# ╔═╡ d8012a6f-b509-44c9-a2b3-eba0513ee037
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

# ╔═╡ 12800da5-a23e-4939-85a7-7ba09f19c501
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


# ╔═╡ f4c8ff31-b33f-4f51-b9d8-8cf15fd8e9b1
KeplerianOrbit(
	rho_s = 0.34,
	R_s = 1.1,
	period = 100.0,
	ecc = 0.3,
	t_0 = 0.5,
	omega = 0.5,
	Omega = 1.0,
	#M_p = 0.1,
	#b = 0.2,
	incl = 0.1,
	M_p = 0.001,
)

# ╔═╡ f01eeef5-24f7-4a31-ab19-a785b78b55e5
orbit = make_orbit_ρₛ(3.0, 0.15^2, 0.3, 0.1)

# ╔═╡ f119776f-9b8c-40bd-a411-a29df4c81a8c
KO(nt::NamedTuple{(:rho_s, :aR_s, :R_s, :period, :ecc, :t_0, :incl, :b, :Omega, :omega, :M_p)}) = 2


# ╔═╡ 33c52bc2-70cb-42f4-87c4-dafd0d3f63fb
@kwcall KO(rho_s=nothing, aR_s=nothing, R_s, period, ecc, t_0, incl=nothing, b=nothing, Omega, omega, M_p=nothing)

# ╔═╡ 91dbdfae-f232-4e95-ab82-9cd9d5bc05f5
KO(
	aR_s = 7.5,
	period = 2.0,
	incl = 90.0,
	t_0 = 0.0,
	ecc = 0.0,
	Omega = 0.0,
	omega = 0.0,
	R_s = 2,
)

# ╔═╡ 3515a6e4-784f-4e1b-baf7-9bd1f4530132
_star_position(orbit, orbit.R_s, 1:10)

# ╔═╡ f4cb6bd7-ae0d-4357-af32-be14c1aceeca
make_orbit_aRₛ() = KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0,
	incl = 90.0,
	t₀ = 0.0,
	ecc = 0.0,
	Ω = 0.0,
	ω = 0.0,
)

# ╔═╡ 5467cc29-d5fb-4d78-9a59-b0f866df9a73
Pl.with_terminal() do 
	@btime make_orbit_aRₛ()
end

# ╔═╡ 07989195-4c8c-4599-ace6-84ad526db2bd
make_orbit_aRₛ()

# ╔═╡ 72506d79-5550-4e10-9830-ccf224bb60bd
make_orbit_ρₛ_units() = KeplerianOrbit(
	ρₛ = 2.0u"g/cm^3",
	Rₛ = 0.5u"Rsun",
	period = 2.0u"d",
	ecc = 0.0,
	t₀ = 0.0u"d",
	incl = 90u"°",
	Ω = 0u"°",
	ω = 0u"°",
	#M_p = 0.5u"Msun"
)

# ╔═╡ 7f40987b-ef60-4786-8639-2100628c4595
make_orbit_aRₛ_units() = KeplerianOrbit(
	aRₛ = 7.5,
	P = 2.0u"d",
	incl = 90.0u"°",
	t₀ = 0.0u"d",
	ecc = 0.0,
	Ω = 0.0u"°",
	ω = 0.0u"°",
)

# ╔═╡ cdda90a8-6f94-4c84-9557-27db687d80e9
# with_terminal() do
# 	@btime make_orbit_ρₛ()
# 	@btime make_orbit_aRₛ()
# end

# ╔═╡ b9a8b205-3df3-4866-bfcf-8841550f9a85
let
	orbit = make_orbit_ρₛ(3.0, 0.15^2, 0.3, 0.1)
	u = [0.4, 0.26] # quad limb dark
	ld = PolynomialLimbDark(u)

	t = range(-1.0, 1.0, length=1000) # days from t0
	rs = range(0, 0.2, length=10) # radius ratio

	fluxes = @. ld(orbit, t, rs')
	
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	for flux in eachcol(fluxes)
		lines!(ax, t, flux)
	end
	
	fig
end

# ╔═╡ a8d30635-7d16-49bf-bf9b-2c6c1e7f8b70
const G_nom = Orbits.G_nom# * u"Rsun^3/Msun/d^2"

# ╔═╡ 1e49a299-a7cc-438a-b186-d2c14619fbf6
ρₛ(P, δ, T, τ) = (3.0*P / (G_nom*π^2)) * (δ^(1/4) / (√(T*τ)))^3

# ╔═╡ 9fe20dfe-6029-45b2-9b48-570f00200ebe
begin
	P = 3.0
	δ = 8.41e-5
	T = 0.5
	τ = 0.125
	
	ρₛ(P, δ, T, τ)
end

# ╔═╡ 39383b74-4c65-4553-afd5-b0d111a9bb03
1.0* u"Rearth/Rsun" |> upreferred

# ╔═╡ 947373bb-c9f6-422c-a200-c9d3e17b49be
let
	ld = IntegratedLimbDark([0.4, 0.26])
	orbit = make_orbit_ρₛ(3.0, 0.15^2, 0.3, 0.1)
	t = range(-1.0, 1.0, length=1000)
	texp = [0.1 0.2 0.3]
	# no extra calculations made
	fluxes = @. ld(orbit, t, 0.2)
	# use quadrature to find time-averaged flux for each t
	fluxes_int = @. ld(orbit, t, 0.2, texp)
	
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	for flux in eachcol(fluxes)
		lines!(ax, t, flux)
	end
	
	for flux in eachcol(fluxes_int)
		lines!(ax, t, flux)
	end
	
	fig
end

# ╔═╡ 2d06d3dd-e266-4da2-a6f5-f5b196189524
let
	ld = SecondaryLimbDark([0.4, 0.26], brightness_ratio=0.1)
	ld_int = IntegratedLimbDark(ld) # composition works flawlessly

	orbit = make_orbit_ρₛ(4.0, 0.15^2, 0.3, 0.1)
	t = range(-1.25, 5, length=1000)
	rs = range(0.01, 0.1, length=6)

	fluxes = @. ld(orbit, t, rs')
	fluxes_int = @. ld_int(orbit, t, rs', texp=0.3)
	
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	for flux in eachcol(fluxes)
		lines!(ax, t, flux)
	end
	
	for flux in eachcol(fluxes_int)
		lines!(ax, t, flux)
	end
	
	fig
end

# ╔═╡ Cell order:
# ╠═afb5a050-b80d-4205-8583-6515779c18e3
# ╠═0895d026-491f-4168-8ad0-509a3227ca89
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═09d609e7-54c4-468c-8f8f-39aa808d20be
# ╠═3de4e3b5-3d64-4e50-b951-de03a8d77ad1
# ╠═72174e0c-0116-4d1e-8365-c1aa63b648ca
# ╠═6de491a0-525d-4283-999f-dff478185214
# ╠═df13cf06-bb11-452f-978e-98d27e226bb2
# ╠═81d987d3-b0e8-48cb-9a9c-d924cd8f7297
# ╠═ca156253-c53d-4b5c-a8b8-7fe694b8f095
# ╠═9d2d78b1-2319-43d2-8a51-c383ba17b585
# ╠═d8012a6f-b509-44c9-a2b3-eba0513ee037
# ╠═12800da5-a23e-4939-85a7-7ba09f19c501
# ╠═2a67f029-8472-47f4-8fdc-054214861355
# ╠═f4c8ff31-b33f-4f51-b9d8-8cf15fd8e9b1
# ╠═5467cc29-d5fb-4d78-9a59-b0f866df9a73
# ╠═f01eeef5-24f7-4a31-ab19-a785b78b55e5
# ╠═f119776f-9b8c-40bd-a411-a29df4c81a8c
# ╠═33c52bc2-70cb-42f4-87c4-dafd0d3f63fb
# ╠═91dbdfae-f232-4e95-ab82-9cd9d5bc05f5
# ╠═3515a6e4-784f-4e1b-baf7-9bd1f4530132
# ╠═f4cb6bd7-ae0d-4357-af32-be14c1aceeca
# ╠═07989195-4c8c-4599-ace6-84ad526db2bd
# ╠═72506d79-5550-4e10-9830-ccf224bb60bd
# ╠═7f40987b-ef60-4786-8639-2100628c4595
# ╠═cdda90a8-6f94-4c84-9557-27db687d80e9
# ╠═b9a8b205-3df3-4866-bfcf-8841550f9a85
# ╠═a8d30635-7d16-49bf-bf9b-2c6c1e7f8b70
# ╠═1e49a299-a7cc-438a-b186-d2c14619fbf6
# ╠═9fe20dfe-6029-45b2-9b48-570f00200ebe
# ╠═39383b74-4c65-4553-afd5-b0d111a9bb03
# ╠═947373bb-c9f6-422c-a200-c9d3e17b49be
# ╠═2d06d3dd-e266-4da2-a6f5-f5b196189524
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
