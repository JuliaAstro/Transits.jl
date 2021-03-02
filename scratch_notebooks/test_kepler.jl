### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 12d523ae-7b24-11eb-3e11-5581a47e1b90
begin
	using Revise
	import Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using PlutoUI
	using StatsPlots
	using Transits.Orbits: relative_position
	using Unitful
	using UnitfulAstro
	using PyCall
end

# ╔═╡ e34d610a-7b23-11eb-0132-55ce740f8e17
pyimport("pip")["main"](["install", "batman-package", "numpy==1.20.1"])

# ╔═╡ e8cf7b42-7b23-11eb-12b9-9978504654a8
begin
	py"""
	import numpy as np
	from batman import _rsky
	
	sky_test = {}
	
	def sky_coords():
		t = np.linspace(-100, 100, 1_000)

		t0, period, a, e, omega, incl = (
			x.flatten()
			for x in np.meshgrid(
				np.linspace(-5.0, 5.0, 2),
				np.exp(np.linspace(np.log(5.0), np.log(50.0), 3)),
				np.linspace(50.0, 100.0, 2),
				np.linspace(0.0, 0.9, 5),
				np.linspace(-np.pi, np.pi, 3),
				np.arccos(np.linspace(0, 1, 5)[:-1]),
			)
		)
		
		r_batman = np.empty((len(t), len(t0)))

		for i in range(len(t0)):
			r_batman[:, i] = _rsky._rsky(
				t, t0[i], period[i], a[i], incl[i], e[i], omega[i], 1, 1
			)
		
		m = r_batman < 100.0
		
		sky_test["m_sum"] = m.sum()
		sky_test["r_batman_m"] = r_batman[m]
		sky_test["m"] = m
		sky_test["t"] = t
		sky_test["t0"] = t0
		sky_test["period"] = period
		sky_test["a"] = a
		sky_test["e"] = e
		sky_test["omega"] = omega
		sky_test["incl"] = incl
	
		return sky_test
	
	"""
	sky_test = py"sky_coords"()
end

# ╔═╡ ec86643a-7b23-11eb-0e9f-7df7963f0452
function compute_r(orbit, t)
	pos = relative_position.(orbit, t)
	r = map(pos) do arr
		√(arr[1]^2 + arr[2]^2)
	end
	return r
end

# ╔═╡ 0be84e52-7b24-11eb-1368-2d1128d1564a
orbits = [
	KeplerianOrbit(
		aRₛ = sky_test["a"][i],
		P = sky_test["period"][i],
		t₀ = sky_test["t0"][i],
		ecc = sky_test["e"][i],
		ω = sky_test["omega"][i],
		incl = sky_test["incl"][i],
	)
	for i in 1:length(sky_test["t0"]) 
]

# ╔═╡ fad76b9a-7b23-11eb-325a-bf2bc3f035c6
begin
	r = Matrix{Float64}(undef, length(sky_test["t"]), length(sky_test["t0"]))
	
	for (orbit, r_i) in zip(orbits, eachcol(r))
		r_i .= compute_r(orbit, sky_test["t"])
	end
	
	r
end

# ╔═╡ b793af16-7b24-11eb-0948-41b16c1e02c6
sum(sky_test["m"]) > 0

# ╔═╡ bcb46be6-7b24-11eb-0de7-0d3d9daaf665
isapprox(r[sky_test["m"]], sky_test["r_batman_m"], atol=2e-5)

# ╔═╡ e757a692-7b24-11eb-2596-c770d1f60fb0
r[sky_test["m"]]

# ╔═╡ ec14e6c2-7b24-11eb-267b-a3a77d1cbd71
sky_test["r_batman_m"]

# ╔═╡ Cell order:
# ╠═12d523ae-7b24-11eb-3e11-5581a47e1b90
# ╠═e34d610a-7b23-11eb-0132-55ce740f8e17
# ╠═e8cf7b42-7b23-11eb-12b9-9978504654a8
# ╠═ec86643a-7b23-11eb-0e9f-7df7963f0452
# ╠═0be84e52-7b24-11eb-1368-2d1128d1564a
# ╠═fad76b9a-7b23-11eb-325a-bf2bc3f035c6
# ╠═b793af16-7b24-11eb-0948-41b16c1e02c6
# ╠═bcb46be6-7b24-11eb-0de7-0d3d9daaf665
# ╠═e757a692-7b24-11eb-2596-c770d1f60fb0
# ╠═ec14e6c2-7b24-11eb-267b-a3a77d1cbd71
