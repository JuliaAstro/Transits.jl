
using LinearAlgebra: dot

struct LimbDarkLightCurve{CT<:AbstractVector,DT}
    u::CT
    c::CT
    c_norm::CT
    driver::DT
end

@doc raw"""
    LimbDarkLightCurve(u)

Quadratically limb-darkened light curve with limb-darkening coefficients `u`. The limb-darkening profile has the quadratic form

```math
\frac{I(\mu)}{I(1)}=1-u_1(1-\mu)-u_2(1-\mu)^2
```
"""
function LimbDarkLightCurve(u::AbstractVector{T}) where T
    u_ext = pushfirst!(u, -one(T))
    c = get_cl(u_ext)
    c_norm = c ./ (π * (c[1] + 2 * c[2] / 3))
    ld = GreensLimbDark(length(c) - 1)
    return LimbDarkLightCurve(u, c, c_norm, ld)
end


function (lc::LimbDarkLightCurve)(orbit::AbstractOrbit, t, r; texp=0, oversample=7, order=0)

    coords = Orbits.relative_position(orbit, t)
    b = sqrt(coords[1]^2 + coords[2]^2)
    if coords[3] > 0
        b_ = abs(b) / orbit.r_star
        r_ = abs(r) / orbit.r_star
        if b_ < 1 + r_
            sT = lc.driver(b_, r_)
            return dot(sT, lc.c_norm) - 1
        end
    end
    return zero(b)
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

struct GreensLimbDark
    lmax::Int
    M
    N

    tmpM
    tmpN
    tmpsT
end

GreensLimbDark(lmax; kwargs...) = GreensLimbDark(Float64, lmax; kwargs...)

function GreensLimbDark(T, lmax::Int; maxiter=100)
    M = greens_M_coeffs(T, lmax; maxiter=maxiter)
    N = greens_N_coeffs(T, lmax; maxiter=maxiter)

    tmpM = similar(M, lmax + 1)
    tmpN = similar(N, lmax + 1)
    tmpsT = similar(M, lmax + 1)
    GreensLimbDark(lmax, M, N, tmpM, tmpN, tmpsT)
end

function (ld::GreensLimbDark)(b, r)
    # complete occultation
    if b < r -1
        @. ld.tmpsT = 0
        return ld.tmpsT
    # no occultation
    elseif iszero(r) || b > r + 1
        @. ld.tmpsT = 0
        ld.tmpsT[1] = π
        ld.tmpsT[2] = 2π / 3
        return ld.tmpsT
    end
    # compute the kite area
    p0, p1, p2 = 1, b, r
    if p0 < p1
        p0, p1 = p1, p0
    end
    if p1 < p2
        p1, p2 = p2, p1
    end
    if p0 < p1
        p0, p1 = p1, p0
    end
    sqarea = (p0 + (p1 + p2)) * (p2 - (p0 - p1)) *
             (p2 + (p0 - p1)) * (p0 + (p1 - p2))
    kite_area2 = sqrt(max(0, sqarea))

    bmr = b - r
    bpr = b + r
    fourbr = 4 * b * r
    onembmr2 = (1 + bmr) * (1 - bmr)
    onembpr2 = (1 + bpr) * (1 - bpr)

    if iszero(b) || iszero(r)
        ksq = Inf
        k = Inf
        kc = 1
        kcsq = 1
        kkc = Inf
        invksq = 0

        ld.tmpsT[1] = π * (1 - r^2)
    else
        ksq = onembpr2 / fourbr + 1
        k = sqrt(ksq)
        if ksq > 1
            kcsq = onembpr2 / onembmr2
            kc = sqrt(kcsq)
            kkc = k * kc
            ld.tmpsT[1] = π * (1 - r^2)
        else
            kcsq = -onembpr2 / fourbr
            kc = sqrt(kcsq)
            kkc = kite_area2 / fourbr
            kap0 = atan(kite_area2, (r - 1) * (r + 1) + b^2)
            kap1 = atan(kite_area2, (1 - r) * (1 + r) + b^2)
            Alens = kap1 + r^2 * kap0 - 0.5 * kite_area2
            ld.tmpsT[1] = π  - Alens
        end
    end

    iszero(ld.lmax) && return ld.tmpsT

    ld.tmpsT[2], Eofk, Em1mKdm = greens_S1(b, r; ksq=ksq, kc=kc)

    ld.lmax == 1 && return ld.tmpsT

    if iszero(b)
        term = 1 - r^2
        dtermdr = -2 * r
        fac = sqrt(term)
        for n in 2:ld.lmax
            ld.tmpsT[n + 1] = -term * r^2 * 2 * π
            term *= fac
        end
        return ld.tmpsT
    end

    # compute the quadratic term
    r2pb2 = (r^2 + b^2)
    η2 = 0.5 * r^2 * (r2pb2 + b^2)
    if ksq > 1
        four_pi_eta = 4 * pi * (η2 - 0.5)
    else
        four_pi_eta = 2 * (-(π - kap1) + 2 * η2 * kap0 - 0.25 * kite_area2 * (1 + 5 * r^2 + b^2))
    end
    ld.tmpsT[3] = 2 * ld.tmpsT[1] * four_pi_eta

    ld.lmax == 2 && return ld.tmpsT

    # now onto the higher order terms
    if ksq < 0.5 && ld.lmax > 3
        downwoardM!(ld.tmpM, ld.M; b=b, r=r, kap0=kap0, ksq=ksq, Eofk=Eofk, Em1mKdm=Em1mKdm, kite_area2=kite_area2)
    else
        upwardM!(ld.tmpM; lmax=ld.lmax, b=b, r=r, kap0=kap0, ksq=ksq, Eofk=Eofk, Em1mKdm=Em1mKdm, kite_area2=kite_area2)
    end

    # computer the remaining terms in the `st` AbstractVector
    ndnp2 = ntuple(n -> n / (n + 2), Val(ld.lmax + 2))
    slice = 4:1 + lmax
    @. ld.tmpsT[slice] = -2 * r^2 * ld.tmpM[slice] + ndnp2[slice] * (
        onemr2mb2 * ld.tmpM[slice] + sqarea * ld.tmpM[slice - 2]
    )

    return ld.tmpsT
end

function greens_S1(b, r; ksq, kc)
    Λ1 = 0
    bmr = b - r
    bpr = b + r
    onembmr2 = (1 - bmr) * (1 + bmr)
    onembpr2 = (1 - bpr) * (1 + bpr)
    if b ≥ 1 + r || iszero(r)
        # no occultation (case 1)
        Λ1 = 0
        Eofk = 0
        Em1mKdm = 0
    elseif b ≤ r - 1
        # full occultation (case 11)
        Λ1 = 0
        Eofk = 0
        Em1mKdm = 0
    else
        if iszero(b)
            # case 10
            sqrt1mr2 = sqrt(1 - r^2)
            Λ1 = -2π * sqrt1mr2^3
            Eofk = 0.5 * π
            Em1mKdm = 0.25 * π
        elseif b ≈ r
            if r ≈ 0.5
                # case 6
                Λ1 = π - 4 / 3
                Eofk = 1
                Em1mKdm = 1
            elseif r < 0.5
                # case 5
                m = 4 * r^2
                Eof = cel(m, 1, 1, 1 - m)
                Em1mKdm = cel(m, 1, 1, 0)
                Λ1 = π + 2 / 3 * ((2 * m - 3) * Eofk - m * Em1mKdm)
            else
                # case 7
                m = r * r^2
                minv = inv(m)
                Eofk = cel(minv, 1, 1, 1 - minv)
                Em1mKdm = cel(minv, 1, 1, 0)
                Λ1 = π + inv(r) / 3 * (-m * Eofk + (2 * m - 3) * Em1mKdm)
            end
        else
            if ksq < 1
                # case 2, case 8
                sqbrinv = 1 / sqbr
                Piofk, Eofk, Em1mKdm = cel(ksq, kc, (b - r)^2 * kc^2, 0, 1, 1, 3 * kc^2 * (b - r) * (b + r), kc^2, 0)
                Λ1 = (1 - mbr) * (1 + mbr) * (Piofk + (-3 + 6 * r^2 + 2 * b * r) * Em1mKdm - 4 * b * r * Eof) * sqbrinv / 3
            elseif ksq > 1
                # case 3, case 9
                bmrdbpr = (b - r) / (b + r)
                μ = 3 * bmrdbpr / onembmr2
                p = bmrdbpr^2 * onembpr2 / onembmr2
                Piofk, Eofk, Em1mKdm = cel(1/ksq, kc, p, 1 + μ, 1, 1, p + μ, kc^2, 0)
                Λ1 = 2 * sqrt(onembmr2) * (onembpr2 * Piofk - (4 - 7 * r^2 - b^2) * Eofk) / 3
            else
                # case 4
                rootr1mr = sqrt(r * (1 - r))
                Λ1 = 2 * acos(1 - 2*r) - 4 / 3 * (3 + 2 * r - 8 * r^2) * rootr1mr - 2 * π * (r > 0.5)
                Eofk = 1
                Em1mKdm = 1
            end
        end
    end

    return ((1 - (r > b)) * 2 * π - Λ1) / 3, Eofk, Em1mKdm
end

function upwardM!(arr; lmax, b, r, sqarea, ksq, kap0, Eofk, Em1mKdm, kite_area2)
    onemr2mb2 = (1 - r) * (1 + r) - b^2

    computeM0123!(arr; b=b, r=r, kap0=kap0, ksq=ksq, Eofk=Eofk, Em1mKdm=Em1mKdm, kite_area2=kite_area2)

    for n in 4:lmax
        arr[n + 1] = (2 * (n - 1) * onemr2mb2 * arr[n - 1] + (n - 2) * sqarea * arr[n - 2]) / n
    end
    return arr
end

function downwardM!(arr, M; b, r, lmax, ksq, sqarea, kap0, Eofk, Em1mKdm, kite_area2, maxiter=100)
    bmr = b - r
    bpr = b + r
    onembmr2 = (1 - bmr) * (1 + bmr)
    sqonembmr2 = sqrt(onembmr2)
    onemr2mb2 = (1 - r) * (1 + r) - b^2
    if ksq < 1
        tol = eps(ksq) * ksq
        term = 0
        fac = k * sqonembmr2^max(0, lmax - 4)
    end
    # now, compute higher order terms until precision reached
    for j in 0:3
        # add leading term to m
        val = M[j + 1, 1]
        k2n = 1

        # compute higher order terms
        for n in 1:maxiter - 1
            k2n *= ksq
            term = k2n * M[j + 1, n + 1]
            val += term
            abs(term) < tol && break
        end
        arr[lmax - 2 + j] = val * fac
        fac *= sqonembmr2
    end

    # recurse downward
    for n in lmax - 4:-1:4
        arr[n + 1] = ((n + r) * arr[n + 5] - 2 * (n + 3) * 
                      onemr2mb2 * arr[n + 3]) / sqarea / (n + 2)
    end

    # compute lowest four exactly    
    computeM0123!(arr; b=b, r=r, kap0=kap0, ksq=ksq, Eofk=Eofk, Em1mKdm=Em1mKdm, kite_area2=kite_area2)

    return arr
end

function computeM0123!(arr; b, r, kap0, ksq, Eofk, Em1mKdm, kite_area2)
    onemr2mb2 = (1 - r) * (1 + r) - b^2
    bmr = b - r
    sqonembmr2 = sqrt((1 + bmr) * (1 - bmr))
    if ksq < 1
        arr[1] = kap0
        arr[2] = 2 * sqbr * 2 * ksq * Em1mKdm
        arr[3] = kap0 * onemr2mb2 + kite_area2
        arr[4] = 8 * sqrt(b * r)^3 * 2 / 3 * ksq * (Eofk + (3 * ksq - 2) * Em1mKdm)
    else
        arr[1] = π
        arr[2] = 2 * sqonembmr2 * Eofk
        arr[3] = π * onemr2mb2
        arr[4] = sqonembmr2^3 * 2 / 3 * ((3 - 2 * invksq) * Eofk + Em1mKdm / ksq)
    end
    return arr
end

greens_M_coeffs(lmax; kwargs...) = greens_M_coeffs(Float64, lmax, kwargs...)
greens_M_coeffs(T, lmax; maxiter=100) = greens_M_coeffs!(Array{T}(undef, 4, maxiter), lmax)
"""
Compute the coefficients of the series expansion for the highest four terms of the `M` integral
"""
function greens_M_coeffs!(out, lmax)
    for j in 0:3
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


greens_N_coeffs(lmax; kwargs...) = greens_N_coeffs(Float64, lmax, kwargs...)
greens_N_coeffs(T, lmax; maxiter=100) = greens_N_coeffs!(Array{T}(undef, 2, maxiter), lmax)

"""
compute the coefficients of the series expansion for the highest two terms of the `N` integral
"""
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
