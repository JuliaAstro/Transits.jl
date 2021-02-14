# Benchmarks

Transits.jl aims to be at least as fast as similar tools. [Limbdark.jl](https://github.com/rodluger/Limbdark.jl) is also written in Julia and [Agol et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..123A/abstract) showed it outperforms starry, PyTransit, and batman in both runtime speed and numerical accuracy. The following benchmarks are works in progress, but they already show a marginal improvement on the Limbdark.jl implementation.

## Setup

!!! warning
    These benchmarks are works in progress

The code can be found in `bench/`. You'll need to set up the environment yourself, including the installation of Limbdark.jl.

## Performance

![](https://github.com/JuliaAstro/Transits.jl/blob/master/bench/timing.png)

## Comparison with Limbdark.jl

![](https://github.com/JuliaAstro/Transits.jl/blob/master/bench/relative_timing.png)

---

![](https://github.com/JuliaAstro/Transits.jl/blob/master/bench/coeff_relative_timing.png)

