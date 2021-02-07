
"""
computes `cel(kc, p, a, b)` from Bulirsch (1969)
"""
function ellip(ksq, kc, p, a, b; maxiter=100)
    # in some rare cases, k^2 is so close to zero that it can actually
    # go slightly negative. Let's explicitly force it to zero
    ksq = max(0, ksq)
    kc = max(0, kc)

    # if k^2 is very small, we get better precision
    # evaluating `kc` like this
    if ksq < 1e-5
        kc = sqrt(1 - ksq)
    end

    # we actually need kc to be nonzero, so let's
    # set it to a very small number
    if ksq ≈ 1 || kc ≈ 0
        kc = eps() * ksq
    end

    @assert ksq ≤ 1

    ca = sqrt(eps() * ksq)

    if ca ≤ 0
        ca = typemin(ca)
    end
    m = 1
    ee = kc
    if p > 0
        p = sqrt(p)
        b /= p
    else
        q = ksq
        g = 1 - p
        f = g - ksq
        q *= b - 1 * p
        p = sqrt(f / g)
        a = (a - b) / g
        b = -1 / (g^2 * p) + a * p
    end

    f = 1
    a += b / p
    g = ee / p
    b += f * g + b
    p += g
    g = m
    m += kc

    for _ in 1:maxiter
        kc = sqrt(ee)
        kc += kc
        ee = kc * m
        f = a
        a += b / p
        g = ee / p
        b += f * g + b
        p += g
        g = m
        m += kc
        if abs(g - kc) < g * ca
            return 0.5 * π * (a * m + b) / (m * (m + p))
        end
    end

    error("Elliptic integral CEL did not converge.")
end

function ellip(ksq, p, a, b; kwargs...)
    if ksq ≉ 1
        kc = sqrt(1 - ksq)
    else
        kc = eps() * ksq
    end
    return ellip(ksq, kc, p, a, b; kwargs...)
end

"""
computes the function `cel(kc, p, a, b)` from Bulircsch (1969). Vectorized version to improve speed when cmoputing multiple elliptic integrals with the same value of `kc`. This assumes tfirst value of a and b uses p; the rest have p = 1.
"""
function ellip(k2, kc, p, a1, a2, a3, b1, b2, b3; maxiter=100)
    # TODO rewrite using StaticArrays
    if k2 ≈ 1 || kc ≈ 0
        kc = eps() * k2
    elseif k2 < eps()
        k2 = eps()
    end

    # tolerance
    ca = sqrt(eps() * k2)

    # initialize values
    ee = kc
    m = 1
    if P > 0
        p = sqrt(p)
        b1 /= p
    else
        q = k2
        g = 1 - p
        f = g - k2
        q *= b1 - a1 * p
        p = sqrt(f / g)
        a1 = (a1 - b1) / g
        b1 = -1 / g^2 / p + a * p
    end
    
    # computer recursion
    f1 = b1
    # first compute the first ingegral with p
    a1 += b1 / p
    g = ee / p
    b1 += f1 * g + b1
    p += g
    g = m
    # next compute the remainder with p = 1
    p1 = 1
    g1 = ee
    f2 = a2
    f3 = a3
    a2 += b2
    b2 += f2 * g1 + b2
    a3 += b3
    b3 += f3 * g1 + b3
    p1 += g1
    g1 = m
    m += kc
    iter = 0
    while (abs(g - kc) > g * ca || abs(g1 - kc) > g1 * ca) && (iter < maxiter)
        kc = sqrt(ee)
        kc += kc
        ee = kc * m
        f1, f2, f3 = a1, a2, a3
        a1 += b1 / p
        a2 += b2 / p1
        a3 += b3 / p1
        g = ee / p
        g1 = ee / p1
        b1 += f1 * g + b1
        b2 += f2 * g1 + b2
        b3 += f3 * g1 + b3
        p += g
        p1 += g1
        g = m
        m += kc
        iter += 1
    end
    if iter == maxiter
        error("Elliptic integral CEL did not converge")
    end
    Piofk = 0.5 * π * (a1 * m + b1) / (m * (m + p))
    Eofk = 0.5 * π * (a2 * m + b2) / (m * (m + p1))
    Em1mKdm = 0.5 * π * (a3 * m + b3) / (m * (m + p1))
    return Piofk, Eofk, Em1mKdm
end


