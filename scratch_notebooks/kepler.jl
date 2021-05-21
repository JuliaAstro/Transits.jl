### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 208b00e9-4d97-4f80-81d1-1b8346b50f83
begin
	using Revise
	using PlutoUI
	using PyCall
	using Transits
	using CairoMakie
end

# ╔═╡ 29f28a74-77f8-11eb-2b70-dd1462a347fc
TableOfContents(depth=6)

# ╔═╡ bddb767e-77f8-11eb-2692-ad86467f0c81
md"## Unitless examples"

# ╔═╡ da4c42b9-a22b-435d-8c1d-c4b93a9cd411
KeplerianOrbit(
		ρₛ = 2.0,
		Rₛ = 0.5,
		P = 2.0,
		ecc = 0.0,
		t₀ = 0.0,
		incl = π / 2.0,
		Ω = 0.0,
		ω = 0.0,
	)

# ╔═╡ dd77f183-4aee-44c3-949b-33ec84747b95
KeplerianOrbit |> methods

# ╔═╡ cc63406d-6154-4052-8788-22ae8dfaaf02
KeplerianOrbit((
	aRₛ = 7.5,
	P = 2.0,
	incl = π / 2.0,
	# b = 0.0,
	t₀ = 0.0,
	ecc = 0.0,
	Ω = 0.0,
	ω = 0.0,
))

# ╔═╡ Cell order:
# ╟─29f28a74-77f8-11eb-2b70-dd1462a347fc
# ╟─bddb767e-77f8-11eb-2692-ad86467f0c81
# ╠═da4c42b9-a22b-435d-8c1d-c4b93a9cd411
# ╠═dd77f183-4aee-44c3-949b-33ec84747b95
# ╠═cc63406d-6154-4052-8788-22ae8dfaaf02
# ╠═208b00e9-4d97-4f80-81d1-1b8346b50f83
