### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using Transits
	using BenchmarkTools
	using PlutoUI
	using CairoMakie
	using Measurements
	using Unitful
	using UnitfulAstro
end

# ╔═╡ 838c7e85-8bd9-4f77-bf95-1f9bf125aced
using JLD2

# ╔═╡ 59419b02-c4b4-4e8c-8870-d5c4924f9245
using Transits.Orbits: relative_position

# ╔═╡ b3e30bfb-695f-44a1-9993-99a0e086f07c
using Test

# ╔═╡ 32b63834-f7c1-4e1c-b6bd-3be0a720dc2a
begin
	stringify_units(value::Unitful.AbstractQuantity, unit) = value
	stringify_units(value, unit) = "$value $unit"
end

# ╔═╡ ce23236f-75a1-4720-8420-f53333227e02
stringify_units(3.0u"Msun/Rsun", "yee")

# ╔═╡ a7dab461-6de1-4cda-8470-2018b3da29e9
0.238732414637843 * (1.0u"Msun/Rsun^3" |> u"g/cm^3").val

# ╔═╡ 87881f39-40c4-49e4-9f4a-6ff56e346742
0.238732414637843u"Msun/Rsun^3" |> u"g/cm^3"

# ╔═╡ a1d0da34-8599-46e8-87e8-47c9512b6140
(1.0u"Msun/Rsun^3" |> u"g/cm^3").val

# ╔═╡ 9acde0d7-cd2e-4e48-b8a0-eea217b09b85
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

# ╔═╡ afb5a050-b80d-4205-8583-6515779c18e3
KeplerianOrbit(
	ρₛ = 2.0u"g/cm^3",
	Rₛ = 1.0u"Rsun",
	P = 2.0u"d",
	ecc = 0.0,
	t₀ = 0.0u"d",
	incl = 90.0u"°",
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

# ╔═╡ 127b92f4-52b0-465f-8679-b93e542a16c0
begin
	sky_coords = load("../test/python_code/test_data/KeplerianOrbit_sky_coords.jld2")
		
	orbits = [
		KeplerianOrbit(
			aRₛ = sky_coords["a"][i],
			P = sky_coords["period"][i],
			incl = sky_coords["incl"][i],
			t₀ = sky_coords["t0"][i],
			ecc = sky_coords["e"][i],
			Ω = 0.0,
			ω = sky_coords["omega"][i],
		)
		for i in 1:length(sky_coords["t0"])
    ]


	# Compute coords
    t = sky_coords["t"]
    x = Matrix{Float64}(undef, length(sky_coords["t"]), length(sky_coords["t0"]))
    y = similar(x)
    z = similar(x)
    for (orbit, x_i, y_i, z_i) in zip(orbits, eachcol(x), eachcol(y), eachcol(z))
        pos = relative_position.(orbit, t) |> as_matrix
        a, b, c = eachcol(pos)
        x_i .= a
        y_i .= b
        z_i .= c
    end
	
	# Compare
    m = sky_coords["m"]
    r = @. √(x^2 + y^2)
    r_Transits = r[m]
    r_batman = sky_coords["r_batman"][m]	
	
# 	fig = Figure(resolution=(1_000, 400))
# 	ax1 = Axis(fig[1, 1])
# 	ax2 = Axis(fig[1, 2])
# # 	#linkyaxes!(ax1, ax2)
	
# 	for (star1, star2, pos) in zip(eachcol(r), eachcol(sky_coords["r_batman"]), ["x", "y"])
# 		lines!(ax1, t, star1, label=pos)
# 		lines!(ax2, t, star2, label=pos)
# 	end
	
# 	axislegend(ax2)
	
# 	fig
end

# ╔═╡ 11452086-26e0-4b4d-96b6-ef7730d027ad
let
i = 1
KeplerianOrbit(
			aRₛ = sky_coords["a"][i],
			P = sky_coords["period"][i],
			incl = sky_coords["incl"][i],
			t₀ = sky_coords["t0"][i],
			ecc = sky_coords["e"][i],
			Ω = 0.0,
			ω = sky_coords["omega"][i],
		)
end

# ╔═╡ 3a6d0079-a02c-4951-a7e5-8b94614a42b9
with_terminal() do
	i = 1
	KeplerianOrbit(
			aRₛ = sky_coords["a"][i],
			P = sky_coords["period"][i],
			incl = sky_coords["incl"][i],
			t₀ = sky_coords["t0"][i],
			ecc = sky_coords["e"][i],
			Ω = 0.0,
			ω = sky_coords["omega"][i],
		)
end

# ╔═╡ 2ce44554-5c01-4beb-a696-6c27b0bfca2d
with_terminal() do
    m_star = 1.3
	r_star = 1.0
	orbit = KeplerianOrbit(
        rho_star = m_star / ((4.0/3.0) * π * r_star^3),
        R_star = r_star,
        period = 100.0,
        ecc = 0.3,
        t_0 = 0.5,
        incl = 0.25*π,
        omega = 0.5,
        Omega = 1.0,
    )
end

# ╔═╡ 4f4ddfec-c9f3-458a-853b-a02c5eadaa8c
zero(0u"m")

# ╔═╡ 81d987d3-b0e8-48cb-9a9c-d924cd8f7297
let
	orbit = KeplerianOrbit(
		M_star = 1.3,
        R_star=1.1,
        t_0 = 0.5,
        period = 100.0,
        ecc = 0.0,
        omega = 0.5,
        Omega = 1.0,
        incl = 0.25 * π,
        M_planet = 0.1,
	)
    
	orbit_flipped = flip(orbit, 0.7)
	
	orbit, orbit_flipped
		
	t = range(0, 100; length=1_000)	

    u_star = _star_position.(orbit, orbit.R_star, t) |> as_matrix
    u_planet_flipped = _planet_position.(orbit_flipped, orbit.R_star, t) |> as_matrix
    for i in 1:3
        @show allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end
	
	println()

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

# ╔═╡ 2d06d3dd-e266-4da2-a6f5-f5b196189524
let
	ld = SecondaryLimbDark([0.4, 0.26], brightness_ratio=0.1)
	ld_int = IntegratedLimbDark(ld) # composition works flawlessly

	orbit = KeplerianOrbit(P = 3.456, ecc=0.0, ω=0.0, incl=π/2, t₀=0.0)
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

# ╔═╡ 9ff128be-c834-45aa-aace-e9e771b390f2


# ╔═╡ daad2939-52b7-455f-a9b0-f5a7a5c52b62
let
   t = range(0, 100; length=1_000)

   orbit = KeplerianOrbit(
       M_star = 1.3,
       M_planet = 0.1,
       R_star = 1.0,
       P = 100.0,
       t_0 = 0.5,
       incl = 45.0,
       ecc = 0.0,
       omega = 0.5,
       Omega = 1.0
   )
   orbit_flipped = flip(orbit, 0.7)

   u_star = as_matrix(_star_position.(orbit, orbit.R_star, t))
   u_planet_flipped = as_matrix(_planet_position.(orbit_flipped, orbit.R_star, t))
   for i in 1:3
       @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
   end

   u_planet = as_matrix(_planet_position.(orbit, orbit.R_star, t))
   u_star_flipped = as_matrix(_star_position.(orbit_flipped, orbit.R_star, t))
   for i in 1:3
       @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
   end
end

# ╔═╡ 40e24355-7451-4b08-b58f-089de5cf39c2
let
    orbit = KeplerianOrbit(
        M_star = 1.3,
        R_star = 1.1,
		t_0 = 0.5,
        period = 100.0,
        ecc = 0.3,
        incl = 0.25*π,
        omega = 0.5,
        Omega = 1.0,
        M_planet = 0.1,
    )

    orbit_flipped = flip(orbit, 0.7)

    t = range(0, 100; length=1_000)

    u_star = as_matrix(_star_position.(orbit, orbit.R_star, t))
    u_planet_flipped = as_matrix(_planet_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_star[:, i], u_planet_flipped[:, i], atol=1e-5)
    end

    u_planet = as_matrix(_planet_position.(orbit, orbit.R_star, t))
    u_star_flipped = as_matrix(_star_position.(orbit_flipped, orbit.R_star, t))
    for i in 1:3
        @test allclose(u_planet[:, i], u_star_flipped[:, i], atol=1e-5)
    end
end

# ╔═╡ Cell order:
# ╠═838c7e85-8bd9-4f77-bf95-1f9bf125aced
# ╠═59419b02-c4b4-4e8c-8870-d5c4924f9245
# ╠═11452086-26e0-4b4d-96b6-ef7730d027ad
# ╠═32b63834-f7c1-4e1c-b6bd-3be0a720dc2a
# ╠═ce23236f-75a1-4720-8420-f53333227e02
# ╠═a7dab461-6de1-4cda-8470-2018b3da29e9
# ╠═87881f39-40c4-49e4-9f4a-6ff56e346742
# ╠═a1d0da34-8599-46e8-87e8-47c9512b6140
# ╠═3a6d0079-a02c-4951-a7e5-8b94614a42b9
# ╠═127b92f4-52b0-465f-8679-b93e542a16c0
# ╠═9acde0d7-cd2e-4e48-b8a0-eea217b09b85
# ╠═afb5a050-b80d-4205-8583-6515779c18e3
# ╠═0895d026-491f-4168-8ad0-509a3227ca89
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═09d609e7-54c4-468c-8f8f-39aa808d20be
# ╠═3de4e3b5-3d64-4e50-b951-de03a8d77ad1
# ╠═72174e0c-0116-4d1e-8365-c1aa63b648ca
# ╠═6de491a0-525d-4283-999f-dff478185214
# ╠═df13cf06-bb11-452f-978e-98d27e226bb2
# ╠═2ce44554-5c01-4beb-a696-6c27b0bfca2d
# ╠═4f4ddfec-c9f3-458a-853b-a02c5eadaa8c
# ╠═81d987d3-b0e8-48cb-9a9c-d924cd8f7297
# ╠═2d06d3dd-e266-4da2-a6f5-f5b196189524
# ╠═9ff128be-c834-45aa-aace-e9e771b390f2
# ╠═daad2939-52b7-455f-a9b0-f5a7a5c52b62
# ╠═40e24355-7451-4b08-b58f-089de5cf39c2
# ╠═b3e30bfb-695f-44a1-9993-99a0e086f07c
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
