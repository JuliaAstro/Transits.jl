using Random: AbstractRNG
using Bijectors
import Bijectors: logabsdetjac, bijector
using Distributions: MultivariateDistribution, Continuous
import Distributions: _rand!, _logpdf
using StatsFuns

"""
    Kipping13()

A non-informative prior for two-parameter limb-darkening coefficients using *triangular sampling* ([Kipping 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.2152K/)).

# Examples

```jldoctest
julia> using Random; rng = Random.seed!(10);

julia> rand(rng, Kipping13())
2-element Vector{Float64}:
 0.24716310305467298
 0.08836997249298882

julia> rand(rng, Kipping13(), 5)
2×5 Matrix{Float64}:
 0.0664907  0.124817   1.00732    0.806902  0.74165
 0.520411   0.222718  -0.389412  -0.314755  0.0768429
```

# References

> [Kipping (2013)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.2152K/)
>
>   "Efficient, uninformative sampling of limb darkening coefficients for two-parameter laws"
"""
struct Kipping13 <: MultivariateDistribution{Continuous} end

Base.length(::Kipping13) = 2

function _rand!(rng::AbstractRNG, ::Kipping13, x::AbstractVector{T}) where {T}
    q1, q2 = rand(rng, T, 2)
    sqrtq1 = sqrt(q1)
    twoq2 = 2 * q2
    x[begin] = sqrtq1 * twoq2
    x[end] = sqrtq1 * (1 - twoq2)
    return x
end

_logpdf(::Kipping13, x::AbstractArray{T}) where {T} = zero(T)

struct Kipping13Transform <: Bijector{1} end

function (::Kipping13Transform)(x::AbstractVector)
    usum = sum(x)
    q = [usum^2, 0.5 * first(x) / usum]
    return @. log(q) - log(1 - q)
end

function (::Inverse{<:Kipping13Transform})(y::AbstractVector)
    tmp = map(logistic, y)
    sqrtq1 = sqrt(first(tmp))
    twoq2 = 2 * last(tmp)
    tmp[begin] = sqrtq1 * twoq2
    tmp[end] = sqrtq1 * (1 - twoq2)
    return tmp
end

logabsdetjac(::Kipping13Transform, y::Number) = -2 * softplus(-y) - y
logabsdetjac(k::Kipping13Transform, y) = sum((yi) -> logabsdetjac(k, yi), y)
bijector(::Kipping13) = Kipping13Transform()
