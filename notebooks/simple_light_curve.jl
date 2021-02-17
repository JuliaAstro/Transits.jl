### A Pluto.jl notebook ###
# v0.12.20

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

# ╔═╡ 68a4138e-70b7-11eb-1ca7-399ee2f63a0e
begin
	import Pkg
	Pkg.activate("..")
	Pkg.instantiate()
	using PlutoUI, Transits, Unitful, UnitfulRecipes, StatsPlots, LaTeXStrings
end

# ╔═╡ 8e62dabc-70c4-11eb-07b3-d7ef7471240b
md"# `Transits.jl` Tutorial 1"

# ╔═╡ 7468b2d0-70c4-11eb-1463-4f3ef693019f
TableOfContents()

# ╔═╡ 8035d886-70c4-11eb-3212-fd9a6cc8de7e
md"## Setup"

# ╔═╡ b32ef55e-70bc-11eb-0c8e-f9f51a88a134
md"""
In this notebook we will explore basic usage of [`Transits.jl`](https://juliaastro.github.io/Transits.jl/stable) by creating some simple light curves. Let's start by first setting up the notebook environment:
"""

# ╔═╡ 9deb2688-70ba-11eb-0222-a7899e86523b
default(palette=:seaborn_colorblind)

# ╔═╡ b4d183f6-70c4-11eb-33be-cd86c368ae35
md"## Orbit definition"

# ╔═╡ 13e902e0-70bd-11eb-00a2-698092b29b2c
md"""
We will now create a `SimpleOrbit` object, which is a subtype of `AbstractOrbit`. This type accepts the following fields:
"""

# ╔═╡ db44287e-70bd-11eb-31ea-a1130ca7077f
fieldnames(SimpleOrbit)

# ╔═╡ f80ed7c2-70bf-11eb-074a-e94f706e5c3b
md"""
Each field is already defined with default values, so we will go ahead and create our `SimpleOrbit` object with a period of 3 days and transit duration of 1 hour, which can be specified seamlessly with the external `Unitful.jl` package: 
"""

# ╔═╡ a0bcc3c4-70b7-11eb-2a20-77644ebb50c1
orbit = SimpleOrbit(period=3u"d", duration=1u"hr")

# ╔═╡ c63edf1e-70c4-11eb-0b86-b17548325707
md"## Limb darkening definition"

# ╔═╡ 65a7b174-70c4-11eb-12bb-1fbc53ee1764
md"""
We can now use this as input into a limb darkening law, which we define next:
"""

# ╔═╡ e7223d2a-70c7-11eb-24e4-af0ff49d42db
u = [0.4, 0.26] # Quadratic limb darkening

# ╔═╡ 570e5624-70b8-11eb-1767-6911283302f5
ld = PolynomialLimbDark(u)

# ╔═╡ 40820d3a-70c5-11eb-04d1-3f32d0d453fe
md"""
For this example, we are using a Polynomial limb darkening law from [Agol, Luger, Foreman-Mackey (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract). This law is part of the `AbstractLimbDark` supertype, which includes the following laws as well that we will explore in other notebooks:
"""

# ╔═╡ 832cd21e-70c5-11eb-05b5-cd9a588751fc
subtypes(Transits.AbstractLimbDark)

# ╔═╡ d0863d52-70c5-11eb-22f6-2fb13f25751c
md"""
The length of the vector of limb darkening coefficients we passed above is equivalent to the order of the polynomial used.
"""

# ╔═╡ 10065930-70c6-11eb-1fa7-551fd12f8e75
md"## Light curve computation"

# ╔═╡ 1fff2420-70c6-11eb-1bc0-9d35bd83a36a
md"""
With the orbit and limb darkening law defined, we can now compute light curves over time `t` for a range of different planet-to-star ratios `rprs`, which can be controlled with the slider below. We perform the computation of each light curve by directly passing our `AbstractOrbit` object to `ld`:
"""

# ╔═╡ 6514b6ea-1374-4852-a4a7-ab8fa496d5d0
md"""
``R_\text{p}/R_* =`` $(@bind rprs Slider(0:0.05:0.2, show_value=true))
"""

# ╔═╡ 8be23852-70b8-11eb-30d3-b315246efe02
begin
	
	# Initialize plot
	p = plot(
		ylim = (0.95, 1.0),
		xguide = "Time from mid-transit",
		yguide = "Relative flux",
		legendtitle = L"R_\mathrm{p}/R_{*}",
		title = "u = $u",
	)
	
	
	# Range of times to compute light curves over
	t = range(-1, 1, length=1_000)u"hr"
	
	# Show previous models
	rprs_previous = 0.0:0.05:rprs
	fluxes_previous = @. ld(orbit, t, rprs_previous')
	plot!(p, t, fluxes_previous, c=:lightgrey, label=false)
	
	# Show current model
	fluxes = @. ld(orbit, t, rprs)
	plot!(t, fluxes, c=1, label=rprs, legend=:bottomleft)
end

# ╔═╡ a5989098-70c7-11eb-24d3-fb1c85b6f770
md"""
Not too bad!
"""

# ╔═╡ af07537e-70bf-11eb-3c0a-592a55780281
note(text) = Markdown.MD(Markdown.Admonition("note", "Note", [text]))

# ╔═╡ ee993496-70bd-11eb-3374-83209e8150de
note(
md"""
These fields are described in the documentation, which can be readily pulled up in the Pluto notebook by selecting the `Live docs` tab below and then clicking on `SimpleOrbit` in the cell above. The `AbstractOrbit` supertype also inlcludes a `KeplerianOrbit` subtype that we will explore in the next notebook.
"""
)

# ╔═╡ 0d935990-70c7-11eb-227a-274d8776916d
note(
md"""
Technically we are passing this to the `compute` function, which is aliased to `ld` for convenience. A row-vector of planet-to-star ratios can be passed as well, and the required light curve computations will be automatically broadcasted:
	
```julia
rprs = [0.0 0.05 0.1 0.15 0.2]
fluxes = @. ld(orbit, t, rprs)
```
"""
)

# ╔═╡ Cell order:
# ╟─8e62dabc-70c4-11eb-07b3-d7ef7471240b
# ╟─7468b2d0-70c4-11eb-1463-4f3ef693019f
# ╟─8035d886-70c4-11eb-3212-fd9a6cc8de7e
# ╟─b32ef55e-70bc-11eb-0c8e-f9f51a88a134
# ╠═68a4138e-70b7-11eb-1ca7-399ee2f63a0e
# ╠═9deb2688-70ba-11eb-0222-a7899e86523b
# ╟─b4d183f6-70c4-11eb-33be-cd86c368ae35
# ╟─13e902e0-70bd-11eb-00a2-698092b29b2c
# ╠═db44287e-70bd-11eb-31ea-a1130ca7077f
# ╟─ee993496-70bd-11eb-3374-83209e8150de
# ╟─f80ed7c2-70bf-11eb-074a-e94f706e5c3b
# ╠═a0bcc3c4-70b7-11eb-2a20-77644ebb50c1
# ╟─c63edf1e-70c4-11eb-0b86-b17548325707
# ╟─65a7b174-70c4-11eb-12bb-1fbc53ee1764
# ╠═e7223d2a-70c7-11eb-24e4-af0ff49d42db
# ╠═570e5624-70b8-11eb-1767-6911283302f5
# ╟─40820d3a-70c5-11eb-04d1-3f32d0d453fe
# ╠═832cd21e-70c5-11eb-05b5-cd9a588751fc
# ╟─d0863d52-70c5-11eb-22f6-2fb13f25751c
# ╟─10065930-70c6-11eb-1fa7-551fd12f8e75
# ╟─1fff2420-70c6-11eb-1bc0-9d35bd83a36a
# ╟─0d935990-70c7-11eb-227a-274d8776916d
# ╟─6514b6ea-1374-4852-a4a7-ab8fa496d5d0
# ╠═8be23852-70b8-11eb-30d3-b315246efe02
# ╟─a5989098-70c7-11eb-24d3-fb1c85b6f770
# ╟─af07537e-70bf-11eb-3c0a-592a55780281
