

struct IntegratedLimbDark{LD<:AbstractLimbDark,T} <: AbstractLimbDark
    limbdark::LD
    texp::T
    order::Int
end

@doc raw"""
    IntegratedLimbDark(limbdark; N=20)

Computes the time-averaged flux in the middle of an exposure by wrapping a limb darkening law `limbdark` with a simple integration. For each time step `t`, `N` extra points are *super-sampled* from `t-texp/2` to `t+texp/2`and the time-averaged flux is calculated using the trapezoidal rule.

**Mathematical form**
```math
\bar{F}(t) = \frac{1}{\Delta t}\int_{t-\Delta t / 2}^{t+\Delta t / 2}{F(t')dt'}
```
where $F$ is the wrapped limb darkening law and $\Delta t$ is the exposure time.
"""
function IntegratedLimbDark(limbdark; N=20)
    return IntegratedLimbDark(limbdark, texp, order)
end
