
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

        s1 = S1()

        if iszero(b)
            term = 1 - r^2
            dtermdr = -2 * r
            fac = sqrt(term)
            dfacdr = -r / fac
            st = zeros(typeof(b), ld.lmax + 1)
            for n in 2:lmax
                st[n + 1] = -term * r^2 * 2 * π
                term *= fac
            end
            return st
        end

        r2pb2 = (r^2 + b^2)
        η2 = 0.5 * r^2 * (r2pb2 + b^2)
        if ksq > 1
            four_pi_eta = 4 * pi * (η2 - 0.5)
        else
            four_pi_eta = 2 * (-(π - kap1) + 2 * η2 * kap0 - 0.25 * kite_area2 * (1 + 5 * r^2 + b^2))
        end
        st[3] = 2 * st[1] * four_pi_eta
        
        lmax == 2 && return st

        # now onto the higher order terms
        if ksq < 0.5 && lmax > 3
            downwoardM()
        else
            upwardM()
        end

        # computer the remaining terms in the `st` AbstractVector
        N = lmax - 2
        @. st[4:4 + N] = -2 * r^2 * m[4:4 + N] + ndnp2[4:4 + N] * (
            onemr2mb2 * m[4:4 + N] + sqarea * m[2:2 + N]
        )



    end

    return st
end

function greens_S1(b, r)
    Λ1 = 0
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
                Eof = ellp(m, 1, 1, 1 - m)
                Em1mKdm = ellip(m, 1, 1, 0)
                Λ1 = π + 2 / 3 * ((2 * m - 3) * Eofk - m * Em1mKdm)
            else
                # case 7
                m = r * r^2
                minv = inv(m)
                Eofk = ellip(minv, 1, 1, 1 - minv)
                Em1mKdm = ellip(minv, 1, 1, 0)
                Λ1 = π + inv(r) / 3 * (-m * Eofk + (2 * m - 3) * Em1mKdm)
            end
        else
            if ksq < 1
                # case 2, case 8
                sqbrinv = 1 / sqbr
                local Piofk
                ellip(ksq, kc, (b - r)^2 * kcsq, 0, 1, 1, 3 * kcsq * (b - r) * (b + r), kcsq, 0, Piofk, Eofk, Em1mKdm)
                Λ1 = (1 - mbr) * (1 + mbr) * (Piofk + (-3 + 6 * r^2 + 2 * b * r) * Em1mKdm - 4 * b * r * Eof) * sqbrinv / 3
            elseif ksq > 1
                # case 3, case 9 
                bmrdbpr = (b - r) / (b + r)
                μ = 3 * bmrdbpr * onembmr2inv
                local Piofk
                ellip(invksq, kc, p, 1 + μ, 1, 1, p + μ, kcsq, 0, Piofk, Eofk, Em1mKdm)
                Λ1 = 2 * sqonembmr2 * (onembpr2 * Piof - (4 - 7 * r^2 - b^2) * Eofk) / 3
            else
                # case 4
                rootr1mr = sqrt(r * (1 - r))
                Λ1 = 2 * acos(1 - 2*r) - 4 / 3 * (3 + 2 * r - 8 * r^2) * rootr1mr - 2 * π * (r > 0.5)
                Eofk = 1
                Em1mKdm = 1
        end

    st = zeros(typeof(b), ld.lmax + 1)
    st[2] = ((1 - (r > b)) * 2 * π - Λ1) / 3

    return st
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
