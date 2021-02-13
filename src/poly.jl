
using SpecialFunctions:loggamma
using LinearAlgebra:dot

struct PolynomialLimbDark{T,VT <: AbstractVector{T},MT <: AbstractMatrix{T},AT <: AbstractArray{T}} <: AbstractLimbDark
    n_max::Int
    u_n::VT
    g_n::VT
    Mn_coeff::AT
    Nn_coeff::MT
    norm::T
    Mn::VT
    Nn::VT
end

@doc raw"""
    PolynomialLimbDark(u::AbstractVector)

Polynomial limb darkening using analytical integrals. The length of the `u` vector is equivalent to the order of polynomial used; e.g., `[0.2, 0.3]` corresponds to quadratic limb darkening.

# Mathematical form
```math
I(\mu) \propto 1 - u_1(1-\mu) - u_2(1-\mu)^2 - \dots - u_N(1-\mu)^N
```
which is equivalent to the series
```math
I(\mu) \propto -\sum_{i=0}^N{u_i(1-\mu)^i}
```
with the definition $u_0 \equiv -1$.
"""
function PolynomialLimbDark(u::AbstractVector{T}; maxiter=100) where T
    # add constant u_0 term
    n_max = length(u)
    u_n = pushfirst!(copy(u), -one(T))
    # get Green's basis coefficients
    g_n = compute_gn(u_n)
    # compute series expansion for M_n and N_n
    Mn_coeff = compute_Mn_coeffs(T, n_max; maxiter=maxiter)
    Nn_coeff = compute_Nn_coeffs(T, n_max; maxiter=maxiter)

    # calculate flux normalization factor, which only depends on first two terms
    if length(g_n) > 1
        norm = inv(π * (g_n[1] + 2 / 3 * g_n[2]))
    else
        norm = inv(π * g_n[1])
    end

    # pre-allocate temp arrays
    Mn = similar(g_n)
    Nn = similar(g_n)

    return PolynomialLimbDark(n_max, u_n, g_n, Mn_coeff, Nn_coeff, norm, Mn, Nn)
end

function compute(ld::PolynomialLimbDark, b::S, r) where S
    T = float(S)

    ## check for trivial cases
    if b ≥ 1 + r || iszero(r)
        # completely unobscured
        return one(T)
    elseif r ≥ 1 + b
        # completely obscured
        return zero(T)
    end

    r2 = r^2
    b2 = b^2

    if iszero(b)
        onemr2 = 1 - r2
        sqrt1mr2 = sqrt(onemr2)

        # annular ellipse
        if ld.n_max == 0
            flux = ld.g_n[begin] * onemr2
        else
            flux = ld.g_n[begin] * onemr2 + 2 / 3 * ld.g_n[begin + 1] * sqrt1mr2^3
        end
        fac = 2 * r2 * onemr2
        for i in 2:ld.n_max
            flux -= ld.g_n[i + 1] * fac
            fac *= sqrt1mr2
        end
        return flux * π * ld.norm
    end

    # take a moment to calculate values used repeatedly in the
    # numerical routines ahead
    onembmr2 = (r + 1 - b) * (1 - r + b)
    onembmr2inv = inv(onembmr2)
    sqonembmr2 = sqrt(onembmr2)
    br = b * r
    fourbr = 4 * br
    fourbrinv = inv(fourbr)
    sqbr = sqrt(br)
    sqbrinv = inv(sqbr)
    onembpr2 = (1 - r - b) * (1 + b + r)
    sqarea = sqarea_triangle(one(T), r, b)
    k2 = max(0, onembpr2 * fourbrinv + 1)
    k = sqrt(k2)
    onemr2mb2 = (1 - r) * (1 + r) - b2
    onemr2pb2 = (1 - r) * (1 + r) + b2
    if k2 > 1
        if k2 > 2
            kc2 = 1 - inv(k2)
        else
            kc2 = onembpr2 * onembmr2inv
        end
    else
        if k2 > 0.5
            kc2 = (r - 1 + b) * (b + r + 1) * fourbrinv
        else
            kc2 = 1 - k2
        end
    end
    kc = sqrt(kc2)

    ## Compute uniform term
    s0, kap0, kite_area2, kck = compute_uniform(b, r; b2, r2, sqarea, fourbrinv)

    flux = ld.g_n[begin] * s0

    if ld.n_max == 0
        return flux * ld.norm
    end

    ## compute linear term
    s1, Eofk, Em1mKdm = compute_linear(b, r; k2, kc, kc2, r2, b2, br, fourbr, sqbr, onembmr2, onembpr2, onembmr2inv)

    flux += ld.g_n[begin + 1] * s1

    if ld.n_max == 1
        return flux * ld.norm
    end

    ## calculate quadratic term
    s2 = compute_quadratic(b, r; s0, r2, b2, kap0, kite_area2, k2)
    flux += ld.g_n[begin + 2] * s2
    if ld.n_max == 2
        return flux * ld.norm
    end

    ## higher order (N > 3) terms require solving Mn and Nn integrals
    if k2 < 0.5 && ld.n_max > 3
        downwardM!(ld.Mn, ld.Mn_coeff; ld.n_max, sqonembmr2, sqbr, onemr2mb2, k, k2, sqarea, kap0, Eofk, Em1mKdm, kite_area2)
    else
        upwardM!(ld.Mn; sqbr, ld.n_max, sqonembmr2, onemr2mb2, sqarea, k2, kap0, Eofk, Em1mKdm, kite_area2)
    end

    # compute remaining terms
    for n in 3:ld.n_max
        pofgn_M = 2 * r2 * ld.Mn[begin + n] - n / (n + 2) *
                       (onemr2mb2 * ld.Mn[begin + n] + sqarea * ld.Mn[begin + n - 2])
        flux -= ld.g_n[begin + n] * pofgn_M
    end

    return flux * ld.norm
end

###
###   Initialization helpers
###

"""
Transform the `u_n` coefficients to `g_n`, which are coefficients in Green's basis from REF.
"""
compute_gn(u::AbstractVector) = compute_gn!(zero(u), u)

function compute_gn!(g_n::AbstractVector{T}, u_n) where T
    n = length(g_n) - 1
    ## First: calculate the a_n terms (eqn. 10)
    a_n = zero(g_n)
    # constant term
    a_n[begin] = one(T)
    for i in 1:n
        bcoeff = one(T) # binomial coefficient
        for j in 0:i
            a_n[begin + j] -= (-1)^j * bcoeff * u_n[begin + i]
            bcoeff *= (i - j) / (j + 1)
        end
    end

    ## Second: calculate the g_n terms (eqn. 15)
    for j in n:-1:2
        if j ≥ n - 1
            g_n[begin + j] = a_n[begin + j] / (j + 2)
        else
            g_n[begin + j] = a_n[begin + j] / (j + 2) + g_n[begin + j + 2]
        end
    end
    if n ≥ 3
        g_n[begin + 1] = a_n[begin + 1] + 3 * g_n[begin + 3]
    elseif n ≥ 1
        g_n[begin + 1] = a_n[begin + 1]
    end
    if n ≥ 2
        g_n[begin] = a_n[begin] + 2 * g_n[begin + 2]
    else
        g_n[begin] = a_n[begin]
    end

    return g_n
end

compute_Mn_coeffs(n_max; kwargs...) = compute_Mn_coeffs(Float64, n_max, kwargs...)
compute_Mn_coeffs(T, n_max; maxiter=100) = compute_Mn_coeffs!(zeros(T, 2, 4, maxiter), n_max)

"""
Compute the series expandsion coefficients for the `M_n` integrals. Need to compute from `n_max - 3:n_max`. Also split for the cases `k^2 < 1` and `k^2 ≥ 1`.
"""
function compute_Mn_coeffs!(Mn_coeff::AbstractArray{T,3}, n_max) where T
    for k2 in 0:1
        coeff = zero(T)
        # loop over n_max to n_max - 3
        for j in 0:3
            # m is n_max to n_max - 3
            m = n_max + j - 3
            mhalf = 0.5 * m
            if k2 < 1
                coeff = sqrt(π) * exp(loggamma(mhalf + 1) - loggamma(mhalf + 1.5))
                Mn_coeff[begin, begin + j, begin] = coeff
                # loop over higher order terms until precision is reached
                for i in 1:size(Mn_coeff, 2) - 1
                    coeff *= (2 * i - 1)^2 / (2 * i * (1 + m + 2 * i))
                    Mn_coeff[begin, begin + j, begin + i] = coeff
                end
            else
                coeff = π
                Mn_coeff[2, begin + j, 1] = coeff
                # loop over higher order terms
                jmax = iseven(m) ? m ÷ 2 : size(Mn_coeff, 2) - 1
                for i in 1:jmax
                    coeff *= (2 + m - 2 * i) * (1 - 2 * i) / (4 * i^2)
                    Mn_coeff[2, begin + j, i + 1] = coeff
                end
            end
        end
    end
    return Mn_coeff
end


compute_Nn_coeffs(n_max; kwargs...) = compute_Nn_coeffs(Float64, n_max, kwargs...)
compute_Nn_coeffs(T, n_max; maxiter=100) = compute_Nn_coeffs!(zeros(T, 2, maxiter), n_max)

"""
Compute the series expandsion coefficients for the `N_n` integrals. Need to compute from `n_max - 3:n_max`. Only need the case `k^2 < 1`.
"""
function compute_Nn_coeffs!(Nn_coeff::AbstractMatrix{T}, n_max) where T
    coeff = zero(T)
    for j in 0:1
        m = n_max + j - 1
        mhalf = 0.5 * m
        coeff = 0.5 * sqrt(π) * exp(loggamma(mhalf + 1) - loggamma(mhalf + 2.5))
        Nn_coeff[begin + j, begin] = coeff
        # now compute higher order terms until precision reached
        for i in 1:size(Nn_coeff, 2) - 1
            coeff *= (4 * i^2 - 1) / (2 * i * (3 + m + 2 * i))
            Nn_coeff[begin + j, begin + i] = coeff
        end
    end

    return Nn_coeff
end

###
###   Computation helpers
###

function sort3(p0, p1, p2)
    if p0 < p1
        p0, p1 = p1, p0
    end
    if p1 < p2
        p1, p2 = p2, p1
    end
    if p0 < p1
        p0, p1 = p1, p0
    end
    return p0, p1, p2
end

function sqarea_triangle(p0, p1, p2)
    p0, p1, p2 = sort3(p0, p1, p2)
    sqarea = (p0 + (p1 + p2)) * (p2 - (p0 - p1)) *
             (p2 + (p0 - p1)) * (p0 + (p1 - p2))
    return sqarea
end

function upwardM!(arr; sqbr, n_max, sqonembmr2, onemr2mb2, sqarea, k2, kap0, Eofk, Em1mKdm, kite_area2)

    Mn_four!(arr; kap0, sqbr, k2, Em1mKdm, onemr2mb2, kite_area2, Eofk, sqonembmr2)

    for n in 4:n_max
        arr[begin + n] = (2 * (n - 1) * onemr2mb2 * arr[begin + n - 2] + (n - 2) * sqarea * arr[begin + n - 4]) / n
    end
    return arr
end

function downwardM!(arr::AbstractVector{T}, Mn_coeff; n_max, sqbr, sqonembmr2, onemr2mb2, k, k2, sqarea, kap0, Eofk, Em1mKdm, kite_area2) where T

    Mn_series!(arr, Mn_coeff; n_max, sqonembmr2, k, k2)

    invsqarea = inv(sqarea)
    # recurse downward
    for n in n_max - 4:-1:4
        arr[begin + n] = ((n + 4) * arr[begin + n + 4] - 2 * (n + 3) *
                      onemr2mb2 * arr[begin + n + 2]) * invsqarea / (n + 2)
    end

    # compute lowest four exactly
    Mn_four!(arr; kap0, sqbr, k2, Em1mKdm, onemr2mb2, kite_area2, Eofk, sqonembmr2)

    return arr
end

function Mn_four!(arr::AbstractVector{T}; kap0, sqbr, k2, Em1mKdm, onemr2mb2, kite_area2, Eofk, sqonembmr2) where T
    # compute lowest four exactly
    if k2 < 1
        arr[begin] = kap0
        arr[begin + 1] = 4 * sqbr * k2 * Em1mKdm
        arr[begin + 2] = kap0 * onemr2mb2 + kite_area2
        arr[begin + 3] = 8 * sqbr^3 * 2/3 * k2 * (Eofk + (3 * k2 - 2) * Em1mKdm)
    else
        arr[begin] = π
        arr[begin + 1] = 2 * sqonembmr2 * Eofk
        arr[begin + 2] = π * onemr2mb2
        arr[begin + 3] = sqonembmr2^3 * 2/3 * ((3 - 2 / k2) * Eofk + Em1mKdm / k2)
    end
end

function Mn_series!(Mn::AbstractVector{T}, Mn_coeff; n_max, sqonembmr2, k, k2) where T
    # Use series expansion to compute M_n:
    # Computing leading coefficient (n=0):
    if k2 < 1
        tol = eps(k2)
        term = zero(T)
        fac = k * sqonembmr2^(n_max - 3)
        # now, compute higher order terms until precision reached
        @inbounds for j in axes(Mn_coeff, 2)
            # add leading term to m
            val = Mn_coeff[begin, j, begin]
            k2n = one(T)
            # compute higher order terms
            for coeff in @view Mn_coeff[begin, j, begin + 1:end]
                k2n *= k2
                term = k2n * coeff
                val += term
                abs(term) < tol && break
            end
            Mn[begin + n_max - 4 + j] = val * fac
            fac *= sqonembmr2
        end
    else # k^2 >= 1
        k2inv = inv(k2)
        tol = eps(k2inv)
        fac = sqonembmr2^(n_max - 3)
        @inbounds for j in axes(Mn_coeff, 2)
            val = Mn_coeff[begin + 1, j, begin]
            k2n = 1
            for coeff in @view Mn_coeff[begin + 1, j, begin + 1:end]
                k2n *= k2inv
                term = k2n * coeff
                val += term
                abs(term) < tol && break
            end
            Mn[begin + n_max - 4 + j] = val * fac
            fac *= sqonembmr2
        end
    end
    return Mn
end

function compute_uniform(b::T, r; sqarea, r2=r^2, b2=b^2, fourbrinv=inv(4*b*r)) where T
    if b ≤ 1 - r
        flux = π * (1 - r2)
        kap0 = convert(T, π)
        kck = zero(T)
        kite_area2 = zero(T)
    else
        kite_area2 = sqrt(sqarea)
        r2m1 = (r - 1) * (r + 1)
        kap0 = atan(kite_area2, r2m1 + b2)
        Πmkap1 = atan(kite_area2, r2m1 - b2)
        kck = kite_area2 * fourbrinv
        flux = Πmkap1 - r2 * kap0 + 0.5 * kite_area2
    end

    return flux, kap0, kite_area2, kck
end

function compute_linear(b::T, r; k2, kc, kc2, onembmr2, onembpr2, onembmr2inv, r2=r^2, b2=b^2, br=b*r, fourbr=4*br, sqbr=sqrt(br)) where T
    if iszero(b) # case 10
        Λ1 = -2 * π * sqrt(1 - r2)^3
        Eofk = 0.5 * π
        Em1mKdm = 0.25 * π
    elseif b == r
        if r == 0.5 # case 6
            Λ1 = π - 4 / 3 - 2 * (b - 0.5) + 6 * (r - 0.5)
            Eofk = one(T)
            Em1mKdm = one(T)
        elseif r < 0.5 # case 5
            m = 4 * r2
            Eofk = cel(m, one(T), one(T), 1 - m)
            Em1mKdm = cel(m, one(T), one(T), zero(T))
            Λ1 = π + 2 / 3 * ((2 * m - 3) * Eofk - m * Em1mKdm) +
                 (b - r) * 4 * r * (Eofk - 2 * Em1mKdm)
        else # case 7
            m = 4 * r2
            minv = inv(m)
            Eofk = cel(minv, one(T), one(T), 1 - minv)
            Em1mKdm = cel(minv, one(T), one(T), zero(T))
            Λ1 = π + 1 / 3 * ((2 * m - 3) * Em1mKdm - m * Eofk) / r -
                 (b - r) * 2 * (2 * Eofk - Em1mKdm)
        end
    else
        if b + r > 1 # case 2, case 8
            Πofk, Eofk, Em1mKdm = cel(k2, kc, (b - r)^2 * kc2, zero(T), one(T), one(T), 3 * kc2 * (b - r) * (b + r), kc2, zero(T))
            Λ1 = onembmr2 * (Πofk + (-3 + 6 * r2 + 2 * br) * Em1mKdm - fourbr * Eofk) / (3 * sqbr)
        elseif b + r < 1 # case 3 case 9
            bmrdbpr = (b - r) / (b + r)
            μ = 3 * bmrdbpr * onembmr2inv
            p = bmrdbpr^2 * onembpr2 * onembmr2inv
            Πofk, Eofk, Em1mKdm = cel(inv(k2), kc, p, 1 + μ, one(T), one(T), p + μ, kc2, zero(T))
            Λ1 = 2 * sqrt(onembmr2) * (onembpr2 * Πofk - (4 - 7 * r2 - b2) * Eofk) / 3
        else # case
            Λ1 = 2 * acos(1 - 2 * r) - 2 * π * (r > 0.5) - (4/3 * (3 + 2 * r - 8 * r2) + 8 * (r + b - 1) * r) * sqrt(r * (1 - r))
            Eofk = one(T)
            Em1mKdm = one(T)
        end
    end

    flux = ((1 - T(r > b)) * 2π - Λ1) / 3

    return flux, Eofk, Em1mKdm
end

function compute_quadratic(b, r;s0, r2, b2, kap0, kite_area2, k2)
    η2 = r2 * (r2 + 2 * b2)
    if k2 ≥ 1
        four_pi_eta = 2 * π * (η2 - 1)
    else
        Πmkap1 = atan(kite_area2, (r - 1) * (r + 1) - b2)
        four_pi_eta = 2 * (-Πmkap1 + η2 * kap0 - 0.25 * kite_area2 * (1 + 5 * r2 + b2))
    end
    return 2 * s0 + four_pi_eta
end