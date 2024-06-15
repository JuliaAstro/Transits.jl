# Benchmarks

Transits.jl aims to be at least as fast as similar tools. [Limbdark.jl](https://github.com/rodluger/Limbdark.jl) is also written in Julia and [ALFM2020](@citet) showed it outperforms starry, PyTransit, and batman in both runtime speed and numerical accuracy. The following benchmarks are works in progress, but they already show a marginal improvement on the Limbdark.jl implementation.

## Setup

!!! warning
    These benchmarks are works in progress

The code can be found in `bench/`. You'll need to set up the environment yourself, including the installation of Limbdark.jl.

## Performance

![](https://raw.githubusercontent.com/JuliaAstro/Transits.jl/main/bench/timing.png)

## Comparison with Limbdark.jl

![](https://raw.githubusercontent.com/JuliaAstro/Transits.jl/main/bench/relative_timing.png)

---

![](https://raw.githubusercontent.com/JuliaAstro/Transits.jl/main/bench/coeff_relative_timing.png)

## References

```@bibliography
Pages = ["bench.md"]
```
