```@meta
CurrentModule = Transits
```

# Transits.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/Transits.jl)
[![Build Status](https://github.com/juliaastro/Transits.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/Transits.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/T/Transits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/Transits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/Transits.jl)

[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4544095.svg)](https://doi.org/10.5281/zenodo.4544095)

Transits.jl provides flexible and powerful occultation curves with limb darkening. The goals of this package are, in this order
* have a simple interface with high *composability*
* be flexible with respect to numeric types and application
* be fully compatible with [ChainRules.jl](https://github.com/juliadiff/ChainRules.jl) automatic differentiation (AD) system to leverage the derived analytical gradients
* provide a codebase that is well-organized, instructive, and easy to extend
* maintain high performance: at least as fast as similar tools

In particular, [`PolynomialLimbDark`](@ref) implements the "starry" limb darkening method, which solves the flux integral analytically. This provides floating-point errors and runtimes that are best in class.

## Installation

To install use [Pkg](https://julialang.github.io/Pkg.jl/v1/managing-packages/). From the REPL, press `]` to enter Pkg-mode

```julia
pkg> add Transits
```
If you want to use the most up-to-date version of the code, check it out from `master`

```julia
pkg> add Transits#master
```

## Citations

If you use Transits.jl or a derivative of it in your work please consider citing it with the [Zenodo DOI](https://doi.org/10.5281/zenodo.4544095). If you use [`PolynomialLimbDark`](@ref) or [`QuadLimbDark`](@ref) please also cite [Agol et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract) and [Luger et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...64L/abstract).
