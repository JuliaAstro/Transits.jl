### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 8c1cd01e-48cd-11ec-1764-9ff6e7f2268f
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using Revise, Transits
end

# ╔═╡ f5e28dec-6358-489f-a5dd-5f19175bba45
using LinearAlgebra

# ╔═╡ e64c66d3-e511-45a6-bef9-69981924997e
orbit = KeplerianOrbit(
	R_star = 0.189,
	M_star = 0.151,
	P = 0.4626413,
	t_0 = 0.2,
	b = 0.5,
	ecc = 0.8,
	omega = 0.1,
)

# ╔═╡ f590bd4c-f89c-46de-b927-0ebca1bb7cac
pos = Transits.Orbits.relative_position.(orbit, orbit.t0)

# ╔═╡ a1d258a3-09e3-4981-a5d6-9c5f515e22f1
hypot(pos[1], pos[2])

# ╔═╡ 2a986329-2b35-4cda-98ec-e4d7ac03edee
norm(pos[1:3])

# ╔═╡ Cell order:
# ╠═8c1cd01e-48cd-11ec-1764-9ff6e7f2268f
# ╠═e64c66d3-e511-45a6-bef9-69981924997e
# ╠═f590bd4c-f89c-46de-b927-0ebca1bb7cac
# ╠═a1d258a3-09e3-4981-a5d6-9c5f515e22f1
# ╠═f5e28dec-6358-489f-a5dd-5f19175bba45
# ╠═2a986329-2b35-4cda-98ec-e4d7ac03edee
