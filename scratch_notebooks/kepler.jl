### A Pluto.jl notebook ###
# v0.14.2

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

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless examples"

# ╔═╡ 22a91084-78a8-11eb-1602-83e6d7cee1bf
make_orbit_ρₛ() = KeplerianOrbit(
	2.0, # ρₛ
	0.5, # Rₛ
	2.0, # P
	0.0, # ecc
	0.0, # t₀
	90.0 * π / 180.0, # incl
)

# ╔═╡ e54167b0-55f3-426a-b7ac-a93e8565297f
orbit_ρₛ = make_orbit_ρₛ()

# ╔═╡ 3e0c5df0-b344-44cc-be6f-b94831a98472
make_orbit_aRₛ() = KeplerianOrbit(
	7.5, # aRₛ
	2.0, # P
	0.0, # b
	0.0, # t₀
	0.0, # ecc
)

# ╔═╡ a6af8efb-ae5e-4f5c-af72-a8687d9bac29
orbit_aRₛ = make_orbit_aRₛ()

# ╔═╡ b7eb5e36-b4a2-43a8-934d-b76adb794f9f
with_terminal() do
	@btime make_orbit_ρₛ()
	@btime make_orbit_aRₛ()
end

# ╔═╡ 8476e7fb-8635-4b91-b52b-2cbf0e1a9ad9
md"""
## Relative positions
"""

# ╔═╡ cdc40849-46b7-4b4a-b78b-bead8c2e7a1d
t = -2:0.05:2 # days

# ╔═╡ 0060649a-24d9-4894-a020-85e850919985
stack(arr_arr) = hcat((reshape(map(p -> p[i], arr_arr), :) for i in 1:3)...)

# ╔═╡ 39bdad37-b58d-40a5-a250-efed63abfc4b
xyz_ρₛ = stack(Transits.Orbits.relative_position.(orbit_ρₛ, t))

# ╔═╡ 8f669a04-0adf-4c0b-90c4-76412aef2332
let
	p = plot(xlabel="Time (d)", ylabel="Positon relative to Rₛ", title="ρₛ")
	plot!(p, t, xyz_ρₛ, label=["x" "y" "z"])
end

# ╔═╡ 984515a8-590e-4e6e-96ca-f1334b8fb6f6
xyz_aRₛ = stack(Transits.Orbits.relative_position.(orbit_aRₛ, t))

# ╔═╡ 237c916f-9434-4c0f-8bf9-3c86d51a12ed
let
	p = plot(xlabel="Time (d)", ylabel="Positon relative to Rₛ", title="aRₛ")
	plot!(p, t, xyz_aRₛ, label=["x" "y" "z"])
end

# ╔═╡ 5274fb9f-999c-440a-aad2-ce5356157b34
md"""
## Light curves
"""

# ╔═╡ f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
u = [0.4, 0.26] # Quadratic LD coeffs

# ╔═╡ 638405b7-0f1e-4eb1-85d6-d17dc04f8699
ld = PolynomialLimbDark(u)

# ╔═╡ 18e62e17-79e0-4586-bfa1-e837c1f17368
fluxes_ρₛ = @. ld(orbit_ρₛ, t, 0.2)

# ╔═╡ 83e060ef-d17d-4561-9808-67a4146642c7
fluxes_aRₛ = @. ld(orbit_aRₛ, t, 0.2)

# ╔═╡ eb9d7c0c-0790-466c-a722-06c9b41b96cd
begin
	p = plot(xlabel="Time (days)", ylabel="Relative flux")
	plot!(p, t, fluxes_ρₛ, lw=5, label="ρₛ")
	plot!(p, t, fluxes_aRₛ, label="aRₛ")
end

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═22a91084-78a8-11eb-1602-83e6d7cee1bf
# ╠═e54167b0-55f3-426a-b7ac-a93e8565297f
# ╠═3e0c5df0-b344-44cc-be6f-b94831a98472
# ╠═a6af8efb-ae5e-4f5c-af72-a8687d9bac29
# ╠═b7eb5e36-b4a2-43a8-934d-b76adb794f9f
# ╟─8476e7fb-8635-4b91-b52b-2cbf0e1a9ad9
# ╠═cdc40849-46b7-4b4a-b78b-bead8c2e7a1d
# ╠═39bdad37-b58d-40a5-a250-efed63abfc4b
# ╠═8f669a04-0adf-4c0b-90c4-76412aef2332
# ╠═984515a8-590e-4e6e-96ca-f1334b8fb6f6
# ╠═237c916f-9434-4c0f-8bf9-3c86d51a12ed
# ╟─0060649a-24d9-4894-a020-85e850919985
# ╟─5274fb9f-999c-440a-aad2-ce5356157b34
# ╠═f8b8bf9b-e667-4c9c-9eb3-126d9b5e71a2
# ╠═638405b7-0f1e-4eb1-85d6-d17dc04f8699
# ╠═18e62e17-79e0-4586-bfa1-e837c1f17368
# ╠═83e060ef-d17d-4561-9808-67a4146642c7
# ╠═eb9d7c0c-0790-466c-a722-06c9b41b96cd
# ╠═fe341026-7278-11eb-2490-0f2ffdeae45e
