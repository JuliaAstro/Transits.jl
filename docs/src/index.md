```@meta
CurrentModule = Transits
```

# Transits.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/Transits/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.org/Transits.jl/dev)

[![CI](https://github.com/JuliaAstro/Transits.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaAstro/Transits.jl/actions/workflows/ci.yml)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/T/Transits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/juliaastro/Transits.jl/graph/badge.svg?token=8kQPWHAWdW)](https://codecov.io/gh/juliaastro/Transits.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4544094.svg)](https://doi.org/10.5281/zenodo.4544094)

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
If you want to use the most up-to-date version of the code, check it out from `main`

```julia
pkg> add Transits#main
```

## Citations

If you use Transits.jl or a derivative of it in your work please consider citing it at the [Zenodo DOI](https://doi.org/10.5281/zenodo.4544094). If you use `PolynomialLimbDark` or `QuadLimbDark` please also cite [Agol et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract) and [Luger et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...64L/abstract). If you use `Kipping13` please cite [Kipping (2013)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.2152K/exportcitation). BibTeX for all those citations can be found in [`CITATIONS.bib`](https://github.com/JuliaAstro/Transits.jl/blob/main/CITATIONS.bib).
