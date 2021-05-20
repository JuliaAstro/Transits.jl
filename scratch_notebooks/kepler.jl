### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

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
	using BenchmarkTools
end

# ╔═╡ 83ed04a9-1913-4b40-afaa-01fa8ecdacd2
using StaticArrays

# ╔═╡ 65dabdc0-b67d-4454-9767-4b9b94d2a11d
using Rotations

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless examples"

# ╔═╡ e54167b0-55f3-426a-b7ac-a93e8565297f
begin
	make_orbit_ρₛ() = KeplerianOrbit((
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	))
	
	make_orbit_ρₛ_KC() = KeplerianOrbit_KC(
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	)
	
	make_orbit_ρₛ_KD() = KeplerianOrbit_KD(
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		#b = 0.0,
		Ω = 0.0,
		ω = 0.0,
	)
	
	orbit_ρₛ,
	orbit_ρₛ_KC,
	orbit_ρₛ_KD = make_orbit_ρₛ(), make_orbit_ρₛ_KC(), make_orbit_ρₛ_KD()
end

# ╔═╡ a6af8efb-ae5e-4f5c-af72-a8687d9bac29
make_orbit_aRₛ() = KeplerianOrbit((
	aRₛ = 7.5,
	P = 2.0,
	incl = π / 2.0,
	# b = 0.0,
	t₀ = 0.0,
	ecc = 0.0,
	Ω = 0.0,
	ω = 0.0,
))

# ╔═╡ f0c2b1b2-b851-4f00-8899-6cf871c190d5
make_orbit_aRₛ()

# ╔═╡ 497322b7-5ff9-4699-8d98-677db5a14061
orbit_aRₛ = make_orbit_aRₛ()

# ╔═╡ b7eb5e36-b4a2-43a8-934d-b76adb794f9f
with_terminal() do
	@btime make_orbit_ρₛ()
	@btime make_orbit_ρₛ_KC()
	@btime make_orbit_ρₛ_KD()
end

# ╔═╡ 5274fb9f-999c-440a-aad2-ce5356157b34
md"""
## Light curves
"""

# ╔═╡ 1a94a407-7a09-4c5f-8796-2570d41ddb4d
t = -0.1:0.01:0.1 # days

# ╔═╡ f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
u = [0.4, 0.26] # Quadratic LD coeffs

# ╔═╡ 638405b7-0f1e-4eb1-85d6-d17dc04f8699
ld = PolynomialLimbDark(u)

# ╔═╡ eb9d7c0c-0790-466c-a722-06c9b41b96cd
begin
	fluxes_ρₛ = @. ld(orbit_ρₛ, t, 0.2)
	fluxes_aRₛ = @. ld(orbit_aRₛ, t, 0.2)
	
	p = plot(xlabel="Time (days)", ylabel="Relative flux")
	plot!(p, t, fluxes_ρₛ, lw=5, label="ρₛ")
	#plot!(p, t, fluxes_aRₛ, label="aRₛ")
end

# ╔═╡ be821fee-6b13-4f02-b7e8-04ed7d56ad39
rotations = [27.014000741644374, 36.63731157445627, 43.95397869463099, 48.503353012894834, 49.99901106593417, 48.34678791380008, 43.650705658060005, 36.20642432789478, 26.482627458951885]

# ╔═╡ 96dc6baf-cc69-4797-a4ac-784feb6897bb
manual = [27.014000741644377, 36.63731157445628, 43.95397869463099, 48.503353012894834, 49.99901106593417, 48.34678791380008, 43.650705658060005, 36.20642432789477, 26.48262745895188]

# ╔═╡ b4dc2635-5a2d-4834-9756-d7a26c007568
separations = 3

# ╔═╡ 6a0e0647-eeb7-491e-b54d-ea8026bad8ed
time = 5:10

# ╔═╡ e5eb611d-e038-40a7-b8d7-b5c957733ce5
function _position(orbit, separation, t)
    sin_ν, cos_ν = Orbits.compute_true_anomaly(orbit, t)
    if iszero(orbit.ecc)
        r = separation
    else
        r = separation * (1 - orbit.ecc^2) / (1 + orbit.ecc * cos_ν)
    end
    X = SA[r * cos_ν, r * sin_ν, zero(r)]
    R = RotZXZ(orbit.ω, -orbit.incl, orbit.Ω)
    return R * X
end

# ╔═╡ 71b217e6-1b19-4e70-b0cd-1ff0637892db
as_matrix(pos) = reinterpret(reshape, Float64, pos) |> permutedims

# ╔═╡ c66c0150-8da0-4ad0-b589-ddc7d1c81a0d
pos_manual = Orbits._position.(orbit_ρₛ, separations, time) |> as_matrix

# ╔═╡ b2b5d53b-8e6d-46c8-93b7-c0e0c9d8a425
pos_rotations = _position.(orbit_ρₛ, separations, time) |> as_matrix

# ╔═╡ 4e1ae817-3d86-4ea2-a66b-1f2b18f67d9b
pos_manual == pos_rotations

# ╔═╡ d9e844b7-d4fa-4ad9-a24f-43da86550560
plot(pos_manual[:, 3] - pos_rotations[:, 3], label="Δz")

# ╔═╡ 1c0a49d5-739c-43b5-b8ee-0ac8d20d4048
plot(rotations - manual)

# ╔═╡ 08b0d3c2-1cee-4a21-ac68-80116d832d9d
function rotate_vector_manual(i, ω, Ω, x, y)
    sin_incl, cos_incl = sincos(i)
    sin_Ω, cos_Ω = sincos(Ω)
    sin_ω, cos_ω = sincos(ω)

    # Rotate about z0 axis by ω
    x1 = cos_ω * x - sin_ω * y
    y1 = sin_ω * x + cos_ω * y

    # Rotate about x1 axis by -incl
    x2 = x1
    y2 = cos_incl * y1
    Z = -sin_incl * y1

    # Rotate about z2 axis by Ω
    X = cos_Ω * x2 - sin_Ω * y2
    Y = sin_Ω * x2 + cos_Ω * y2

    return SA[X, Y, Z]
end

# ╔═╡ ad1d9770-8eeb-4b3c-a8e2-2c73be3c9b74
function rotate_vector_Rotations(i, ω, Ω, x, y)
	X = SA[x, y, zero(x)]
    R = RotZXZ(Ω, -i, ω)
	return R * X
end

# ╔═╡ 20f770bd-b9ed-4d8e-a7df-f31b24226e1a
inputs = π/3, 0, 0.1, 3.7, 4

# ╔═╡ 04412ecf-d080-4fbf-aabf-4727743b311a
rotate_vector_manual(inputs...)

# ╔═╡ ca1474c6-66e6-494c-9a69-971036fa8099
rotate_vector_Rotations(inputs...)

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═e54167b0-55f3-426a-b7ac-a93e8565297f
# ╠═a6af8efb-ae5e-4f5c-af72-a8687d9bac29
# ╠═f0c2b1b2-b851-4f00-8899-6cf871c190d5
# ╠═497322b7-5ff9-4699-8d98-677db5a14061
# ╠═b7eb5e36-b4a2-43a8-934d-b76adb794f9f
# ╟─5274fb9f-999c-440a-aad2-ce5356157b34
# ╠═1a94a407-7a09-4c5f-8796-2570d41ddb4d
# ╠═f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
# ╠═638405b7-0f1e-4eb1-85d6-d17dc04f8699
# ╠═eb9d7c0c-0790-466c-a722-06c9b41b96cd
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
# ╠═be821fee-6b13-4f02-b7e8-04ed7d56ad39
# ╠═96dc6baf-cc69-4797-a4ac-784feb6897bb
# ╠═b4dc2635-5a2d-4834-9756-d7a26c007568
# ╠═6a0e0647-eeb7-491e-b54d-ea8026bad8ed
# ╠═c66c0150-8da0-4ad0-b589-ddc7d1c81a0d
# ╠═b2b5d53b-8e6d-46c8-93b7-c0e0c9d8a425
# ╠═4e1ae817-3d86-4ea2-a66b-1f2b18f67d9b
# ╠═d9e844b7-d4fa-4ad9-a24f-43da86550560
# ╠═e5eb611d-e038-40a7-b8d7-b5c957733ce5
# ╠═83ed04a9-1913-4b40-afaa-01fa8ecdacd2
# ╠═71b217e6-1b19-4e70-b0cd-1ff0637892db
# ╠═1c0a49d5-739c-43b5-b8ee-0ac8d20d4048
# ╠═65dabdc0-b67d-4454-9767-4b9b94d2a11d
# ╠═08b0d3c2-1cee-4a21-ac68-80116d832d9d
# ╠═ad1d9770-8eeb-4b3c-a8e2-2c73be3c9b74
# ╠═20f770bd-b9ed-4d8e-a7df-f31b24226e1a
# ╠═04412ecf-d080-4fbf-aabf-4727743b311a
# ╠═ca1474c6-66e6-494c-9a69-971036fa8099
