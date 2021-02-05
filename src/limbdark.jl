
struct LimbDarkLightCurve{CT<:AbstractVector}
    u::CT
    c::CT
    c_norm::CT
end

function LimbDarkLightCurve(u::AbstractVector{T}) where T
    u_ext = pushfirst!(u, -one(T))
    c = get_cl(u_ext)
    c_norm = c ./ (π * (c[1] + 2 * c[2] / 3))
    return LimbDarkLightCurve(u, c, c_norm)
end


function (lc::LimbDarkLightCurve)(orbit::AbstractOrbit, t, r; texp=0, oversample=7, order=0)

    coords = relative_position(orbit, t)
    b = sqrt(coords[1]^2 + coords[2]^2)
    f = limbdark(lc.c_norm, b, r, coords[3])
end


function get_cl(u::AbstractVector{T}) where T
    N = length(u)

    c = similar(u)

    a = zero(u)
    a[firstindex(u)] = one(T)
    
    # compute the a_n coefficients
    for i in 1:N - 1
        bcoeff = 1
        sign = 1
        for j in 0:i
            a[j + 1] -= u[i + 1] * bcoeff * sign
            sign *= -1
            bcoeff *= (i - j) / (j + 1)
        end
    end

    # now compute the c_n coefficients
    for j in N-1:-1:max(2, N-2)
        c[j + 1] = a[j + 1] / (j + 1)
    end
    
    for j in N-3:-1:2
        c[j + 1] = a[j + 1] / (j + 2) + c[j + 2]
    end

    if N ≥ 4
        c[2] = a[2] + 3 * c[4]
    else
        c[2] = a[2]
    end

    if N ≥ 3
        c[1] = a[1] + 2 * c[3]
    else
        c[1] = a[1]
    end

    return c
end

function limbdark(cl::AbstractVector{T}, b, r, los) where T
    ld = GreensLimbDark(length(cl))
    if los > 0
        b_ = abs(b)
        r_ = abs(r)
        if b_ < 1 + r_
            st = ld(b_, r_)
            return st * cl - 1 
        end
    end
    return zero(T)
end

struct GreensLimbDark
    lmax::Int
    b
    r
    k
    kc
    kkc
    kap0
    kap1
    invksq
    fourbr
    invfourbr
    bmr
    bpr
    onembmr2
    omembmr2inv
    sqonembmr2
    onembpr2
    b2mr22
    onemr2mb2
    sqarea
    sqbr
    kite_area2
    Eofk
    Em1mKdm

    M
    N
    M_coeffs
    N_coeffs

    n_
    invn
    ndnp2
end

function GreensLimbDark(lmax::Int)
    computeMCoeff()
    computeNCoeff()
    for n in 0 : lmax + 2

    end
    GreensLimbDark(lmax, )
end

function (ld::GreensLimbDark)(b, r)
    # complete occultation
    if b < r -1
        return zeros(typeof(b), ld.lmax + 1)
    # no occultation
    elseif iszero(r) || b > r + 1
        st = zeros(typeof(b), ld.lmax + 1)
        st[1] = π
        st[2] = 2π / 3
        return st
    else
        # compute the kite area
        p0, p1, p2 = 1, b, r
        p0 < p1 && (p0, p1 = p1, p0)
        p1 < p2 && (p1, p2 = p2, p1)
        p0 < p1 && (p0, p1 = p1, p0)
        sqarea = (p0 + (p1 + p2)) * (p2 - (p0 - p1)) * (p2 + (p0 - p1)) * (p0 + (p1 - p2))
        kite_area2 = sqrt(max(0, sqarea))

        bmr = b - r
        bpr = b + r

        if iszero(b) || iszero(r)
            ksq = Inf
            k = Inf
            kc = 1

            kcsq = 1

            kkc = Inf
            invksq = 0
            st = zeros(typeof(b), ld.lmax + 1)
            st[1] = π * (1 - r^2)
        else
            bpr = b + r
            ksq = (1 + bpr) * (1 - bpr) / (4 * b * r) + 1
            k = sqrt(ksq)
            if ksq > 1
                kcsq = (1 + bpr) * (1 - bpr) / (1 + bmr) * (1 - bmr)
                kc = sqrt(kcsq)
                kkc = k * kc
                st = zeros(typeof(b), ld.lmax + 1)
                st[1] = π * (1 - r^2)
            else
                kcsq = -(1 + bpr) * (1 - bpr) / (4 * b * r)
                kc = sqrt(kcsq)
                kkc = kite_area2 / (4 * b * r)
                kap0 = atan(kite_area2, (r - 1) * (r + 1) + b^2)
                kap1 = atan(kite_area2, (1 - r) * (1 + r) + b^2)
                Alens = kap1 + r^2 * kap0 - 0.5 * kite_area2
                st = zeros(typeof(b), ld.lmax + 1)
                st[1] = π  - Alens
            end

            end
        end


    end
end

greens_M_coeffs(lmax; kwargs...) = greens_M_coeffs(Float64, lmax, kwargs...)
greens_M_coeffs(T, lmax; maxiter=100) = greens_M_coeffs!(Array{T}(undef, 4, maxiter), lmax)

function greens_M_coeffs!(out, lmax)
    for j in 0:4
        n = lmax - 3 + j
        coeff = sqrt(π) * wallis(n)
        out[j + 1, 1] = coeff
        # now compute higher order terms until precision reached
        for i in 1:size(out, 2) - 1
            coeff *= (2 * i - 1)^2 / (2 * i * (1 + n + 2 * i))
            out[j + 1, i + 1] = coeff
        end
    end
    return out
end

"""
    wallis(n)

Compute the Wallis ratio recursively `n` times.

```math
\\Gamma(1 + n/2) / \\Gamma(3/2 + n/2)
```
"""
function wallis(n)
    if iseven(n)
        z = 1 + n/2
        dz = -1
    else
        z = 1 + (n - 1)/2
        dz = 0
    end
    A = 1.0
    B = sqrt(π)
    for i in 1:z + dz - 1
        A *= i + 1
        B *= i - 0.5
    end
    for i in max(1, z+dz):z
        B *= i - 0.5
    end
    return iseven(n) ?  A / B : B / a
end


greens_N_coeffs(lmax; kwargs...) = greens_N_coeffs(Float64, lmax, kwargs...)
greens_N_coeffs(T, lmax; maxiter=100) = greens_N_coeffs!(Array{T}(undef, 2, maxiter), lmax)

function greens_N_coeffs!(out::AbstractMatrix{T}, lmax) where T
    coeff = zero(T)
    # ksq < 1
    for j in 0:1
        n = lmax - 1 + j

        # add leading term to N
        coeff = sqrt(π) * wallis(n) / (n + 3)
        out[j + 1, 1] = coeff

        # now computer higher order terms until precision is reached
        for i in 1:size(out, 2) - 1
            coeff *= (4 * i^2 - 1) / (2 * i * (3 + n + 2 * i))
            out[j + 1, i + 1] = coeff
        end

    end

    return out
end
