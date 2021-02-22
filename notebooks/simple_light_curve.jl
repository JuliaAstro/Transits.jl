### A Pluto.jl notebook ###
# v0.12.21

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
	default(palette=:seaborn_colorblind)
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

# ╔═╡ a24d5b2a-754e-11eb-2bf5-3d5b7f4a0ec6
md"""
Now we need to gather the ingredients that make up a transit: 1) an orbit, and 2) a limb-darkening law. Let's start with describing the orbit of our exoplanetary system.
"""

# ╔═╡ b4d183f6-70c4-11eb-33be-cd86c368ae35
md"## Orbit definition"

# ╔═╡ 13e902e0-70bd-11eb-00a2-698092b29b2c
md"""
We need to specify the parameters that define our orbital system. These parameters are stored in instances of the `AbstractOrbit` supertype. The most straightforward example of this is the `SimpleOrbit` subtype, so we will focus on using this subtype for the rest of this learning notebook. Let's see how to create this object next.
"""

# ╔═╡ f80ed7c2-70bf-11eb-074a-e94f706e5c3b
md"""
From the live docs, we see that `SimpleOrbit` accepts the following parameters:

* `period` - The orbital period of the planets, nominally in days

* `duration` - The duration of the transit, same units as period

* `t0` - The midpoint time of the reference transit, same units as period

* `b` - The impact parameter of the orbit

* `r_star` - The radius of the star, nominally in solar radii

Likewise, we see that `t0`, `b`, and `r_star` already have default values defined. The `;` in the function signature `SimpleOrbit(; period, duration, t0=0, b=0, r_star=1)` means that everything to the right is a keyword argument that we can define by name, so let's go ahead and do that to define an orbit with the following parameters:

* `period` = 3 days
* `duration` = 1 hour
* `t0` = 0 days
* `b` = 0
* `r_star` = 1 (in units of solar radii)
"""

# ╔═╡ a0bcc3c4-70b7-11eb-2a20-77644ebb50c1
orbit = SimpleOrbit(period=3u"d", duration=1u"hr")

# ╔═╡ 3086d2ae-7554-11eb-3638-2f94c1011bab
md"""
With our orbit defined, we now turn to our second and last ingredient: the limb darkening law.
"""

# ╔═╡ c63edf1e-70c4-11eb-0b86-b17548325707
md"## Limb darkening definition"

# ╔═╡ a9f8fa38-7554-11eb-36ad-a76d087d95c6
md"""
Limb darkening is the natural phenomenon where a star's surface brightness appears to fall off as we look from its center to off-angle towards its limbs. The rate at which it falls off is described by a given limb darkening law. For this example, we are using a polynomial limb darkening law from [Agol, Luger, Foreman-Mackey (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract). This law uses analytically derived integrals, which makes it very fast and numerically accurate. There are other limb darkening laws which share a similar `AbstractLimbDark` interface. This interface lets us create composite laws like `IntegratedLimbDark` or `SecondaryLimbDark` with almost no effort. We will explore these laws and more in further notebooks.

Pulling up the documentation for our particular law, `PolynomialLimbDark`, we see that it just accepts a single parameter `u`, the vector of limb darkening coefficients. These are defined such that:

```math
I(\mu) \propto -\sum_{i=0}^N u_i(1 - \mu)^i, \quad u_0 \equiv -1\quad,
```

where ``μ`` is a dimensionless parameter that varies from ``1`` to ``0`` as we move from the center to the edge of the star, as seen projected on the sky, ``N`` is the order of the polynomial, ``I`` is the intensity of the star, and ``u_i`` are the components of `u`.

For this example, we will use a quadratic limb darkening law ``(N = 2)``:

```math
\begin{align}
I(\mu) &\propto -u_0(1 - \mu)^0 - u_1(1 - μ) - u_2(1 - μ)^2 \\
&= 1 - u_1(1 - μ) - u_2(1 - μ)^2 \quad.
\end{align}
```

"""

# ╔═╡ ba86f45c-755b-11eb-29b3-f7e6163857d9
md"""
With the groundwork set, we can now return to the creation of our `PolynomialLimbDark` object, for a quadratic limb darkening law defined by ``\mathbf u = [0.4, 0.26]``:
"""

# ╔═╡ e7223d2a-70c7-11eb-24e4-af0ff49d42db
u = [0.4, 0.26] # Quadratic limb darkening

# ╔═╡ 570e5624-70b8-11eb-1767-6911283302f5
ld = PolynomialLimbDark(u)

# ╔═╡ 10065930-70c6-11eb-1fa7-551fd12f8e75
md"## Light curve computation"

# ╔═╡ 1fff2420-70c6-11eb-1bc0-9d35bd83a36a
md"""
With the orbit and limb darkening law defined, we can now compute light curves over time `t` for a range of different planet-to-star ratios `rprs`, which can be controlled with the slider below. We perform the computation of each light curve by directly passing our `AbstractOrbit` object to `ld`:
"""

# ╔═╡ b6a523e6-754b-11eb-112c-d3b1787d9818
begin
	rprs_range = 0:0.05:0.2
	cs = Dict()
	for (i ,rprs) in enumerate(rprs_range)
		cs[rprs] = i
	end
end

# ╔═╡ 6514b6ea-1374-4852-a4a7-ab8fa496d5d0
md"""
``R_\text{p}/R_* =`` $(@bind rprs Slider(rprs_range, show_value=true))
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
	plot!(t, fluxes, c=cs[rprs], label=rprs, legend=:bottomleft)
end

# ╔═╡ a5989098-70c7-11eb-24d3-fb1c85b6f770
md"""
Not too bad!
"""

# ╔═╡ af07537e-70bf-11eb-3c0a-592a55780281
tip(text) = Markdown.MD(Markdown.Admonition("tip", "Tip", [text]))

# ╔═╡ ee993496-70bd-11eb-3374-83209e8150de
tip(
md"""
Information about functions/types/etc. can be pulled up directly in Pluto by selecting the `Live docs` tab in the bottom right corner and then clicking on the item that you would like to explore. These live docs can also automatically be called by typing `?` in a cell, followed directly by what you would like to search for.
"""
)

# ╔═╡ 9bb70f6a-7553-11eb-3ba2-5f70e3e63c11
tip(
md"""
Parameters with units can be specified seamlessly with external [`Unitful.jl`](https://painterqubits.github.io/Unitful.jl/stable/) package using the `u"<unit here>"` string macro, and everything should just work™.
"""
)

# ╔═╡ 1e7360b4-755b-11eb-0fff-9f0445a6bbd7
note(text) = Markdown.MD(Markdown.Admonition("note", "Note", [text]))

# ╔═╡ 27746f32-755b-11eb-2c0f-e7eea21bee93
note(md"""
The total number of terms is always one more than the order of the polynomial (i.e., the length of our limb darkening coefficient vector) because of the convention ``u_0 \equiv -1`` that we set. This ensures that the first term in our polynomial limb darkening law is ``1``, which corresponds to the maximum apparent brightness that the star can achieve.
"""
)

# ╔═╡ 0d935990-70c7-11eb-227a-274d8776916d
note(
md"""
All of the limb darkening laws can be called like a function: `ld(b, r)`, or using the `compute` method: `compute(ld, b, r)`. We provide both for convenience. It's natural and easy to treat the laws as functions, but some code is easier to write using the concrete `compute` method. In addition, we don't provide any default vectorized versions of the methods, since Julia does not need to be vectorized for speed (unlike numpy). This lets users apply our laws in a variety of ways (e.g., `map`, `broadcast`, `Threads.@threads`, etc.). Here is an example using [broadcasting](https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting) with the `@.` macro to create a grid of outputs across the `t` and `rprs` column and row vectors:
	
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
# ╟─a24d5b2a-754e-11eb-2bf5-3d5b7f4a0ec6
# ╟─b4d183f6-70c4-11eb-33be-cd86c368ae35
# ╟─13e902e0-70bd-11eb-00a2-698092b29b2c
# ╟─ee993496-70bd-11eb-3374-83209e8150de
# ╟─f80ed7c2-70bf-11eb-074a-e94f706e5c3b
# ╟─9bb70f6a-7553-11eb-3ba2-5f70e3e63c11
# ╠═a0bcc3c4-70b7-11eb-2a20-77644ebb50c1
# ╟─3086d2ae-7554-11eb-3638-2f94c1011bab
# ╟─c63edf1e-70c4-11eb-0b86-b17548325707
# ╟─a9f8fa38-7554-11eb-36ad-a76d087d95c6
# ╟─27746f32-755b-11eb-2c0f-e7eea21bee93
# ╟─ba86f45c-755b-11eb-29b3-f7e6163857d9
# ╠═e7223d2a-70c7-11eb-24e4-af0ff49d42db
# ╠═570e5624-70b8-11eb-1767-6911283302f5
# ╟─10065930-70c6-11eb-1fa7-551fd12f8e75
# ╟─1fff2420-70c6-11eb-1bc0-9d35bd83a36a
# ╟─0d935990-70c7-11eb-227a-274d8776916d
# ╟─6514b6ea-1374-4852-a4a7-ab8fa496d5d0
# ╟─b6a523e6-754b-11eb-112c-d3b1787d9818
# ╠═8be23852-70b8-11eb-30d3-b315246efe02
# ╟─a5989098-70c7-11eb-24d3-fb1c85b6f770
# ╟─af07537e-70bf-11eb-3c0a-592a55780281
# ╟─1e7360b4-755b-11eb-0fff-9f0445a6bbd7
