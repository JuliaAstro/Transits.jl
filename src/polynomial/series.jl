#= 
This file contains the code for calculating the Mn and Nn integrals via
series expansion. 
=#

###
### Mn Series
###

compute_Mn_coeffs(n_max; kwargs...) = compute_Mn_coeffs(Float64, n_max, kwargs...)
compute_Mn_coeffs(T, n_max; maxiter=100) = compute_Mn_coeffs!(zeros(T, 2, 4, maxiter), n_max)

"""
Compute the series expandsion coefficients for the `M_n` integrals. Need to compute from `n_max - 3:n_max`. Also split for the cases `k^2 < 1` and `k^2 ≥ 1`.
"""
function compute_Mn_coeffs!(Mn_coeff::AbstractArray{T,3}, n_max) where T
    for k2 in (0, 1)
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
                for i in 1:size(Mn_coeff, 3) - 1
                    coeff *= (2 * i - 1)^2 / (2 * i * (1 + m + 2 * i))
                    Mn_coeff[begin, begin + j, begin + i] = coeff
                end
            else
                coeff = π
                Mn_coeff[begin + 1, begin + j, begin] = coeff

                # loop over higher order terms
                jmax = iseven(m) ? m ÷ 2 : size(Mn_coeff, 3) - 1
                for i in 1:jmax
                    coeff *= (2 + m - 2 * i) * (1 - 2 * i) / (4 * i^2)
                    Mn_coeff[begin + 1, begin + j, begin + i] = coeff
                end
            end
        end
    end
    return Mn_coeff
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
    for m in n_max - 4:-1:4
        arr[begin + m] = ((m + 4) * arr[begin + m + 4] - 2 * (m + 3) *
                      onemr2mb2 * arr[begin + m + 2]) * invsqarea / (m + 2)
    end

    # compute lowest four exactly
    Mn_four!(arr; kap0, sqbr, k2, Em1mKdm, onemr2mb2, kite_area2, Eofk, sqonembmr2)

    return arr
end

function Mn_four!(arr::AbstractVector{T}; kap0, sqbr, k2, Em1mKdm, onemr2mb2, kite_area2, Eofk, sqonembmr2) where T
    # compute lowest four exactly
    if k2 ≤ 1 # eqn 77,79
        arr[begin] = kap0
        arr[begin + 1] = 4 * sqbr * k2 * Em1mKdm
        arr[begin + 2] = kap0 * onemr2mb2 + kite_area2
        arr[begin + 3] = 8 * sqbr^3 * 2 / 3 * k2 * (Eofk + (3 * k2 - 2) * Em1mKdm)
    else # eqn 78,80
        arr[begin] = π
        arr[begin + 1] = 2 * sqonembmr2 * Eofk
        arr[begin + 2] = π * onemr2mb2
        arr[begin + 3] = sqonembmr2^3 * 2 / 3 * ((3 - 2 / k2) * Eofk + Em1mKdm / k2)
    end
    return arr
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

###
### Nn Series
###

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


function Nn_lower!(Nn::AbstractVector, Nn_coeff; Mn, onembpr2, sqonembmr2, n_max, k, k2, kwargs...)
    Nn_series!(Nn, Nn_coeff; k2, sqonembmr2, n_max, k)

    for m in n_max - 2:-1:2
        Nn[begin + m] = ((m + 4) * Nn[begin + m + 2] - Mn[begin + m + 2]) / (onembpr2 * (m + 2))
    end

    Nn_two!(Nn; k, k2, kwargs...)
    return Nn
end

function Nn_series!(Nn::AbstractVector{T}, Nn_coeff; k2, sqonembmr2, n_max, k) where T
    tol = eps(k2)
    term = zero(T)
    fac = sqonembmr2^(n_max - 1) * k * k2
    for j in 1:2
        Nn[begin + n_max - 2 + j] = Nn_coeff[begin + j - 1, begin]
        k2n = one(T)
        for n in 1:size(Nn_coeff, 2) - 1
            k2n *= k2
            term = k2n * Nn_coeff[begin + j - 1, begin + n]
            Nn[begin + n_max - 2 + j] += term
            abs(term) < tol && break
        end
        Nn[begin + n_max - 2 + j] *= fac
        fac *= sqonembmr2
    end

    return Nn
end

function Nn_raise!(Nn::AbstractVector, Mn; n_max, onembpr2, kwargs...)
    Nn_two!(Nn; kwargs...)
    for m in 2:n_max
        Nn[begin + m] = (Mn[begin + m] + m * onembpr2 * Nn[begin + m - 2]) / (m + 2)
    end
    return Nn
end

function Nn_two!(Nn::AbstractVector; kap0, kc, k2, sqbr, k, Eofk, Em1mKdm)
    if k2 ≤ 1
        Nn[begin] = 0.5 * kap0 - k * kc
        Nn[begin + 1] = 4 / 3 * sqbr * k2 * (2 * Em1mKdm - Eofk)
    else
        Nn[begin] = 0.5 * π
        Nn[begin + 1] = 4 / 3 * sqbr * k * (2 * Eofk - Em1mKdm)
    end

    return Nn
end
