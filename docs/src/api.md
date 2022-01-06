# API/Reference

## Index

```@index
```

### Light Curves

```@docs
AbstractLimbDark
(AbstractLimbDark)(args...; kwargs...)
PolynomialLimbDark
QuadLimbDark
IntegratedLimbDark
SecondaryLimbDark
compute
compute(::AbstractLimbDark, ::Orbits.AbstractOrbit, t, r)
```

### Gradients

Gradients and jacobians are integrated directly into [ChainRules.jl](https://github.com/juliadiff/ChainRules.jl) via `frule`s and `rrule`s. *For most users, this just means using AD libraries like [ForwardDiff.jl](https://github.com/juliadiff/ForwardDiff.jl) and [Zygote.jl](https://github.com/FluxML/Zygote.jl) is effortless and fast*.

```jldoctest grads
using Transits
using Zygote

lightcurve(X) = compute(PolynomialLimbDark(X[3:end]), X[1], X[2])
grad(X) = lightcurve'(X) # Zygote gradient
grad([0.1, 0.1, 0.4, 0.26])

# output
4-element Vector{Float64}:
  0.0004972185834858653
 -0.2419262730830416
 -0.0048107583897073185
 -0.0024501564976671724
```

To help demonstrate the logic behind these chain rules, here we derive a simple gradient function manually.

```jldoctest grads
using ChainRules

u_n = [0.4, 0.26]
μ = 0.1
ror = 0.1
X0 = [μ, ror, u_n...]

function gradr(X)
    ld, ld_pullback = rrule(PolynomialLimbDark, X[3:end])
    f, f_pullback = rrule(compute, ld, X[1], X[2])

    f̄ = one(eltype(X))
    _, l̄d, b̄, r̄ = f_pullback(f̄)
    _, ū_n = ld_pullback(l̄d)
    return [b̄, r̄, ū_n...]
end

gradr([0.1, 0.1, 0.4, 0.26])

# output
4-element Vector{Float64}:
  0.0004972185834858653
 -0.2419262730830416
 -0.0048107583897073185
 -0.0024501564976671724
```

For the most granular support for gradients and jacobians, peer into the depths of `polynomial/poly-grad.jl` and `polynomial/quad-grad.jl`. These functions are not part of the public API and are not guaranteed any stability according to [semantic versioning](https://semver.org/).

### Orbits

```@docs
SimpleOrbit
KeplerianOrbit
Orbits.relative_position
```

### Distributions

```@docs
Kipping13
```
