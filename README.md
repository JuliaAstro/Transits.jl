# Transits.jl

[![Build Status](https://github.com/juliaastro/Transits.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/Transits.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/T/Transits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/Transits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/Transits.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/Transits.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/Transits.jl/dev)

Flexible photometric transit curves with limb darkening

**WIP**: Currently under progress by @mileslucas

## Usage

```julia
using Transits

orbit = SimpleOrbit(period=3, duration=0.5)
u = [0.4, 0.26] # quad limb dark
ld = PolynomialLimbDark(u)

t = range(-1, 1, length=1000) # days from t0
rs = range(0, 0.2, length=10) # radius ratio

fluxes = @. ld(orbit, t, rs')
```

```julia
using ColorSchemes, Plots
plot(t, fluxes, xlabel="time - t0 [d]", ylabel="relative flux",
     leg=false, title="Quadratic Limb Darkening (u=[0.4, 0.26])",
     palette=palette(:inferno, size(fs, 2)*2))
```

![](limbdark.png)
