
using FastGaussQuadrature

struct IntegratedLimbDark{LD<:AbstractLimbDark,WT,NT} <: AbstractLimbDark
    driver::LD
    order::Int
    nodes::NT
    weights::WT
end

@doc raw"""
    IntegratedLimbDark(limbdark; N=21, basis=:legendre)
    IntegratedLimbDark(u; kwargs...)

Computes the time-averaged flux in the middle of an exposure by wrapping a limb darkening law `limbdark` with a quadrature scheme. For each time step `t`, `N` extra points are *super-sampled* from `t-texp/2` to `t+texp/2`and the time-averaged flux is calculated via quadrature.

If a set of limb darkening coefficients, `u`, is provided, a [`PolynomialLimbDark`](@ref) law will be used by default.

# Mathematical form
```math
\bar{F}(t) = \frac{1}{\Delta t}\int_{t-\Delta t / 2}^{t+\Delta t / 2}{F(t')dt'}
```
where $F$ is the wrapped limb darkening law and $\Delta t$ is the exposure time.

## Quadrature
The integration is approximated via [Guassian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature)
```math
\frac{1}{\Delta t} \int{F(t')dt'} \approx \frac12\sum_i^N{w_i * F(\frac{\Delta t}{2}\xi_i + t)}
```
where the weights `w_i` and nodes `ξ_i` are defined by the given quadrature rule. The nodes are defined by evaluating orthogonal polynomials `N` times between -1 and 1. Notice the change of interval required to go from the natural bounds of the orthogonal polynomial basis, `-1, 1`, to the range defined by the exposure time.

The following bases are available from [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl). In addition, a function can be passed which calculates `nodes, weights = f(N)`.
* `:legendre` - Legendre polynomial base on the open `(-1, 1)`
* `:radau` - Legendre polynomial base on the semi-open `[-1, 1)` interval
* `:lobatto` - Legendre polynomial base on the closed `[-1, 1]` interval
"""
function IntegratedLimbDark(driver::AbstractLimbDark; N=21, basis=:legendre)
    nodes, weights = get_nodes_weights(basis, N)
    return IntegratedLimbDark(driver, N, nodes, weights)
end

get_nodes_weights(basis, N) = basis(N)

function get_nodes_weights(basis::Symbol, N)
    if basis === :legendre
        return gausslegendre(N)
    elseif basis === :radau
        return gaussradau(N)
    elseif basis === :lobatto
        return gausslobatto(N)
    else
        error("$basis not a support quadrature basis")
    end
end

function IntegratedLimbDark(u::AbstractVector; kwargs...)
    driver = PolynomialLimbDark(u)
    return IntegratedLimbDark(driver; kwargs...)
end

compute(ld::IntegratedLimbDark, orbit::AbstractOrbit, t, r; texp=nothing) = 
    compute(ld, orbit, t, r, texp)

function compute(ld::IntegratedLimbDark, orbit::AbstractOrbit, t, r, texp)
    # perform change of interval
    half_texp = 0.5 * texp
    # perform quadrature
    flux = sum(zip(ld.weights, ld.nodes)) do (w, ξ)
        f = compute(ld.driver, orbit, half_texp * ξ + t, r)
        return w * f
    end
    return 0.5 * flux
end

compute(ld::IntegratedLimbDark, orbit::AbstractOrbit, t, r, ::Nothing) = 
    compute(ld.driver, orbit, t, r)

function Base.show(io::IO, ld::IntegratedLimbDark{L1}) where L1
    N = length(ld.nodes)
    np = L1.name.name
    print(io, "IntegratedLimbDark($np, $N)")
end

function Base.show(io::IO, ::MIME"text/plain", ld::IntegratedLimbDark)
    p = ld.driver
    N = length(ld.nodes)
    print(io, "IntegratedLimbDark\n driver: $p\n N: $N")
end