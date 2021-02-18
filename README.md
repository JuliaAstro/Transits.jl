# Transits.jl

[![Build Status](https://github.com/juliaastro/Transits.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/Transits.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/T/Transits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/Transits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/Transits.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/Transits.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/Transits.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4544095.svg)](https://doi.org/10.5281/zenodo.4544095)

Flexible photometric transit curves with limb darkening. The goals of this package are, in this order

* have a simple interface with high *compasibility*
* be flexible with respect to numeric types and application
* be fully compatible with [ChainRules.jl](https://github.com/juliadiff/ChainRules.jl) automatic differentiation (AD) system to leverage the derived analytical gradients
* provide a codebase that is well-organized, instructive, and easy to extend
* maintain high performance: at least as fast as similar tools

**WIP**: Currently under progress by @mileslucas

### Current TODO list

in some order of importance

- [ ] Finish writing `KeplerOrbit` (help wanted)
- [ ] Gradients using ChainRulesCore
- [ ] Gradient tests using ChainRulesTestUtils
- [ ] Kipping prior distributions (with gradients) (help wanted)
- [ ] documenation section "Introduction" describing transits and talking about science, very expository (help wanted)
- [ ] Plotting functinos (recreate ALFM 20 plots)
    * recipe for lightcurve which automatically makes a simple orbit and shows features
    * similar as above but special one for SecondaryLimbDark to side-by-side plot secondary
    * figure 3 and 6 can be written with recipe
- [ ] examples (show rich composability of julia)
    * We can use Pluto notebooks for examples that are *learning-oriented*
    * For tutorials and *problem-oriented* examples prefer a Julia script that can be weaved into the docs (with Literate.jl e.g.) (or just as easily weaved into a jupyter notebook!)
- [ ] benchmarks (recreate ALFM 20 plots)
    * I have some code in `bench/speed.jl`. This needs modularized- the `benchmark` function can be rewritten like `benchmark(f, Nu, Npts)` and abstracted for `f` as types.
    * once the code is reorganized (maybe even put the new `benchmark` in a module for future thinking) decide whether to save data or save images, then build the link with the documentation
    * need to add benchmarks for the python code. PyCall has worked fine for benchmarking for me in the past, so lets write the driver in Julia. [here is a link](https://github.com/rodluger/Limbdark.jl/blob/master/tex/figures/python/compare_to_batman.py) to the plotting code from ALFM
    * in general, the more this can be automated the better (including CI!)
- [ ] look at simpson integrated light curve (ALFM 20)


If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/Transits.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/Transits.jl/discussions) and join or open a new topic.

## Installation

To install use [Pkg](https://julialang.github.io/Pkg.jl/v1/managing-packages/). From the REPL, press `]` to enter Pkg-mode

```julia
pkg> add Transits
```
If you want to use the most up-to-date version of the code, check it out from `master`

```julia
pkg> add Transits#master
```

## Usage

```julia
using Transits

orbit = SimpleOrbit(period=3, duration=1)
u = [0.4, 0.26] # quad limb dark
ld = PolynomialLimbDark(u)

t = range(-1, 1, length=1000) # days from t0
rs = range(0, 0.2, length=10) # radius ratio

fluxes = @. ld(orbit, t, rs')
```

![](limbdark.png)

## Integrated and Secondary Curves

`IntegratedLimbDark` can be used to numerically integrate each light curve exposure in time

```julia
ld = IntegratedLimbDark([0.4, 0.26])
orbit = SimpleOrbit(period=3, duration=1)
t = range(-1, 1, length=1000)
texp = [0.1 0.2 0.3]
# no extra calculations made
flux = @. ld(orbit, t, 0.2)
# use quadrature to find time-averaged flux for each t
flux_int = @. ld(orbit, t, 0.2, texp)
```

![](integrated.png)

`SecondaryLimbDark` can be used to generate secondary eclipses given a surface brightness ratio

```julia
ld = SecondaryLimbDark([0.4, 0.26], brightness_ratio=0.1)
ld_int = IntegratedLimbDark(ld) # composition works flawlessly

orbit = SimpleOrbit(period=4, duration=1)
t = range(-1.25, 2.75, length=1000)
rs = range(0.01, 0.1, length=6)

f = @. ld(orbit, t, rs')
f_int = @. ld_int(orbit, t, rs', texp=0.3)
```

![](secondary.png)

## Using Units

Units from `Unitful.jl` are a drop-in substitution for numbers

```julia
using Unitful
orbit = SimpleOrbit(period=10u"d", duration=5u"hr")
t = range(-6, 6, length=1000)u"hr"
flux = @. ld(orbit, t, 0.1)
```

## Citations

If you use Transits.jl or a derivative of it in your work please consider citing it at the [Zenodo DOI](https://doi.org/10.5281/zenodo.4544095). If you use `PolynomialLimbDark` or `QuadLimbDark` please also cite [Agol et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract) and [Luger et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...64L/abstract).
