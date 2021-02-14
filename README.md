# Transits.jl

[![Build Status](https://github.com/juliaastro/Transits.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/Transits.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/T/Transits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/Transits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/Transits.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/Transits.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/Transits.jl/dev)

Flexible photometric transit curves with limb darkening. The goals of this package are, in this order

* have a simple interface with high *compasibility*
* be flexible with respect to numeric types and application
* be fully compatible with [ChainRules.jl](https://github.com/juliadiff/ChainRules.jl) automatic differentiation (AD) system to leverage the derived analytical gradients
* provide a codebase that is well-organized, instructive, and easy to extend
* maintain high performance: at least as fast as similar tools

**WIP**: Currently under progress by @mileslucas

### Current TODO list before v0.1.0

- [x] ~~finish writing unit tests~~
- [x] ~~fix numerical errors in Mn integral~~ (Thanks so much for the help @icweaver)
- [ ] prepare DOI ~~and references for starry/Agol/etc.~~

Current TODOs that are further down the horizon, in some order of importance

- [ ] Finish writing `KeplerOrbit`
- [ ] Gradients using ChainRulesCore
- [ ] Gradient tests using ChainRulesTestUtils
- [ ] Plotting functinos (recreate ALFM 20 plots)
- [ ] examples (show rich composability of julia)
- [ ] benchmarks (recreate ALFM 20 plots)
- [ ] look at simpson integrated light curve (ALFM 20)

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/Transits.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/Transits.jl/discussions) and join or open a new topic.

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

If you use Transits.jl or a derivative of it in your work please consider citing it (DOI to be made). If you use `PolynomialLimbDark` or `QuadLimbDark` please also cite [Agol et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract) and [Luger et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...64L/abstract).
