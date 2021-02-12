```@meta
CurrentModule = Transits
```

# Transits.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/Transits.jl)
[![Build Status](https://github.com/juliaastro/Transits.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/Transits.jl/actions)
[![Coverage](https://codecov.io/gh/juliaastro/Transits.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/Transits.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Transits.jl provides flexible and powerful occultation curves with limb darkening. The goals of this package are, in this order
* have a simple interface with high *compasibility*
* be as flexible with respect to numeric types and application
* be fully compatible with [ChainRules.jl](https://github.com/juliadiff/ChainRules.jl) automatic differentiation (AD) system to leverage the derived analytical gradients
* provide a code-base that is well-organized, instructive, and easy to extend
* maintain high performance: at least as fast as similar tools

In particular, [`PolynomialLimbDark`](@ref) implements the "starry" limb darkening method, which solves the flux integral analytically. This provides floating-point errors and runtimes that are best in class.

