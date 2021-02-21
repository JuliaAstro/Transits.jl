import ChainRulesCore: frule, rrule


function frule((_, _), ::typeof(compute_gn), u_n::AbstractVector{T}) where T
    compute_gn_jac(u_n)
end

function compute_gn_jac(u_n::AbstractVector{T}) where T
    n = length(u_n) - 1
    g_n = zero(u_n)
    dgdu = fill!(similar(u_n, n + 1, n + 1), zero(T))

    ## First: calculate the a_n terms (eqn. 10)
    a_n = zero(g_n)
    dadu = copy(dgdu)
    # constant term
    a_n[begin] = one(T)
    for i in 1:n
        bcoeff = one(T) # binomial coefficient
        for j in 0:i
            r = (-1)^j * bcoeff
            a_n[begin + j] -= r * u_n[begin + i]
            dadu[begin + j, begin + i] -= r
            bcoeff *= (i - j) / (j + 1)
        end
    end

    ## Second: calculate the g_n terms (eqn. 15)
    for j in n:-1:2
        if j ≥ n - 1
            g_n[begin + j] = a_n[begin + j] / (j + 2)
            for i in 0:n - 1
                dgdu[begin + j, begin + i + 1] = dadu[begin + j, begin + i + 1] / (j + 2)
            end
        else
            g_n[begin + j] = a_n[begin + j] / (j + 2) + g_n[begin + j + 2]
            for i in 0:n - 1
                dgdu[begin + j, begin + i + 1] = dadu[begin + j, begin + i + 1] / (j + 2) + dgdu[begin + j + 2, begin + i + 1]
            end
        end
    end
    if n ≥ 3
        g_n[begin + 1] = a_n[begin + 1] + 3 * g_n[begin + 3]
        for i in 0:n - 1
            dgdu[begin + 1, begin + i + 1] = dadu[begin + 1, begin + i + 1] + 3 * dgdu[begin + 3, begin + i + 1]
        end
    elseif n ≥ 1
        g_n[begin + 1] = a_n[begin + 1]
        for i in 0:n - 1
            dgdu[begin + 1, begin + i + 1] = dadu[begin + 1, begin + i + 1]
        end
    end
    if n ≥ 2
        g_n[begin] = a_n[begin] + 2 * g_n[begin + 2]
        for i in 0:n - 1
            dgdu[begin, begin + i + 1] = dadu[begin, begin + i + 1] + 2 * dgdu[begin + 2, begin + i + 1]
        end
    else
        g_n[begin] = a_n[begin]
        for i in 0:n - 1
            dgdu[begin, begin + i + 1] = dadu[begin, begin + i + 1]
        end
    end

    return g_n, dgdu
end

function compute_gn_jac(u_n::StaticVector{3,T}) where T
    g_n = compute_gn(u_n)
    ∇g_n = SA[zero(T) -one(T)  -1.5
              zero(T)  one(T)   2.0
              zero(T)  zero(T) -0.25]
    return g_n, ∇g_n
end

function compute_uniform_grad(b::T, r; r2, b2, sqarea, fourbrinv) where T
    if b ≤ 1 - r
        flux = π * (1 - r2)
        kap0 = convert(T, π)
        kck = zero(T)
        kite_area2 = zero(T)
        ∇s0 = SA[zero(T), -2 * π * r]
    else
        kite_area2 = sqrt(sqarea)
        r2m1 = (r - 1) * (r + 1)
        kap0 = atan(kite_area2, r2m1 + b2)
        Πmkap1 = atan(kite_area2, r2m1 - b2)
        kck = kite_area2 * fourbrinv
        flux = Πmkap1 - r2 * kap0 + 0.5 * kite_area2
        ∇s0 = SA[kite_area2 / b, -2 * r * kap0]
    end

    return (flux, kap0, kite_area2, kck), ∇s0
end

function frule((_, _, _), ::typeof(compute_uniform), b, r; kwargs...)
    compute_uniform_grad(b, r; kwargs...)
end


function compute_linear_grad(b::T, r; k2, kc, kc2, onembmr2, onembpr2, onembmr2inv, sqonembmr2, r2=r^2, b2=b^2, br=b*r, fourbr=4*br, sqbr=sqrt(br)) where T
    if iszero(b) # case 10
        Λ1 = -2 * π * sqrt(1 - r2)^3
        Eofk = 0.5 * π
        Em1mKdm = 0.25 * π

        ∇s1 = SA[zero(T), -2 * π * r * sqrt1mr2]
    elseif b == r
        if r == 0.5 # case 6
            Λ1 = π - 4 / 3 - 2 * (b - 0.5) + 6 * (r - 0.5)
            Eofk = one(T)
            Em1mKdm = one(T)
            ∇s1 = SA[2/3, -2]
        elseif r < 0.5 # case 5
            m = 4 * r2
            Eofk = cel(m, one(T), one(T), 1 - m)
            Em1mKdm = cel(m, one(T), one(T), zero(T))
            Λ1 = π + 2 / 3 * ((2 * m - 3) * Eofk - m * Em1mKdm) +
                 (b - r) * 4 * r * (Eofk - 2 * Em1mKdm)
            ∇s1 = SA[-4/3 * r * (Eofk -2 * Em1mKdm), -4 * r * Eofk]
        else # case 7
            m = 4 * r2
            minv = inv(m)
            Eofk = cel(minv, one(T), one(T), 1 - minv)
            Em1mKdm = cel(minv, one(T), one(T), zero(T))
            Λ1 = π + 1 / 3 * ((2 * m - 3) * Em1mKdm - m * Eofk) / r -
                 (b - r) * 2 * (2 * Eofk - Em1mKdm)
            ∇s1 = SA[2/3 * (2 * Eofk - Em1mKdm), -2 * Em1mKdm]
        end
    else
        if b + r > 1 # case 2, case 8
            Πofk, Eofk, Em1mKdm = cel(k2, kc, (b - r)^2 * kc2, zero(T), one(T), one(T), 3 * kc2 * (b - r) * (b + r), kc2, zero(T))
            Λ1 = onembmr2 * (Πofk + (-3 + 6 * r2 + 2 * br) * Em1mKdm - fourbr * Eofk) / (3 * sqbr)
            ∇s1 = SA[2 * r * onembmr2 * (-Em1mKdm + 2 * Eofk) / (sqbr * 3), -2 * r * onembmr2 * Em1mKdm / sqbr]
        elseif b + r < 1 # case 3 case 9
            bmrdbpr = (b - r) / (b + r)
            μ = 3 * bmrdbpr * onembmr2inv
            p = bmrdbpr^2 * onembpr2 * onembmr2inv
            Πofk, Eofk, Em1mKdm = cel(inv(k2), kc, p, 1 + μ, one(T), one(T), p + μ, kc2, zero(T))
            Λ1 = 2 * sqrt(onembmr2) * (onembpr2 * Πofk - (4 - 7 * r2 - b2) * Eofk) / 3
            ∇s1 = SA[-4/3 * r * sqonembmr2 * (Eofk - 2 * Em1mKdm), -4 * r * sqonembmr2 * Eofk]
        else # case
            Λ1 = 2 * acos(1 - 2 * r) - 2 * π * (r > 0.5) - (4/3 * (3 + 2 * r - 8 * r2) + 8 * (r + b - 1) * r) * sqrt(r * (1 - r))
            Eofk = one(T)
            Em1mKdm = one(T)
            g = -8 * r * sqrt(r * (1 - r))
            ∇s1 = SA[-g/3, g]
        end
    end

    flux = ((1 - T(r > b)) * 2π - Λ1) / 3

    return (flux, Eofk, Em1mKdm), ∇s1
end

function frule((_, _, _), ::typeof(compute_linear), b, r; kwargs...)
    compute_linear_grad(b, r; kwargs...)
end


function compute_quadratic_grad(b, r; s0, r2, b2, kap0, kite_area2, k2, ∇s0)
    r2pb2 = r2 + b2
    η2 = r2 * (r2pb2 + b2)
    if k2 > 1
        four_pi_eta = 2 * π * (η2 - 1)
        ∇η = 4 * π * SA[2 * b * r2, 2 * r * r2pb2]
    else
        Πmkap1 = atan(kite_area2, (r - 1) * (r + 1) - b2)
        four_pi_eta = 2 * (-Πmkap1 + η2 * kap0 - 0.25 * kite_area2 * (1 + 5 * r2 + b2))
        ∇η = SA[2/b * (4 * b2 * r2 * kap0 - (1 + r2pb2) * kite_area2),
                8 * r * (r2pb2 * kap0 - kite_area2)]
    end
    s2 = 2 * s0 + four_pi_eta
    ∇s2 = 2 * ∇s0 + ∇η
    return  s2, ∇s2
end

function frule((_, _, _), ::typeof(compute_quadratic), b, r; kwargs...)
    compute_quadratic_grad(b, r; kwargs...)
end


function compute_grad(ld::PolynomialLimbDark, b::S, r) where S
    T = float(S)
    bcut = 1e-3
    n = ld.n_max
    dfdg = zero(ld.g_n)

    ## check for trivial cases
    if b ≥ 1 + r || iszero(r)
        # completely unobscured
        return one(T), dfdg, @SVector zeros(T, 2)
    elseif r ≥ 1 + b
        # completely obscured
        return zero(T), dfdg, @SVector zeros(T, 2)
    end

    r2 = r^2
    b2 = b^2

    if iszero(b)
        onemr2 = 1 - r2
        sqrt1mr2 = sqrt(onemr2)

        # annular ellipse
        facd = -2 * r
        if ld.n_max == 0
            flux = ld.g_n[begin] * onemr2
            dfdrb[1] = ld.g_n[begin] * facd
        else
            flux = ld.g_n[begin] * onemr2 + 2 / 3 * ld.g_n[begin + 1] * sqrt1mr2^3
            dfdrb[1] = ld.g_n[begin] * facd + ld.g_n[begin + 1] * facd * sqrt1mr2
        end
        fac = 2 * r2 * onemr2
        for i in 2:ld.n_max
            flux -= ld.g_n[begin + i] * fac
            dfdrb[1] += ld.g_n[begin + i] * facd * (1 * onemr2 - i * r2)
            dfdg[begin + i] -= fac
            fac *= sqrt1mr2
            facd *= sqrt1mr2
        end
        dfdrb *= π * ld.norm
        dfdg[1] = (onemr2 - flux) * π * ld.norm
        dfdg[2] = 2/3 * (sqrt1mr2^3 - flux) * π * ld.norm
        return flux * π * ld.norm, dfdg * ld.norm, dfdrb
    end

    # take a moment to calculate values used repeatedly in the
    # numerical routines ahead
    onembmr2 = (r - b + 1) * (1 - r + b)
    onembmr2inv = inv(onembmr2)
    sqonembmr2 = sqrt(onembmr2)
    br = b * r
    fourbr = 4 * br
    fourbrinv = inv(fourbr)
    sqbr = sqrt(br)
    sqbrinv = inv(sqbr)
    onembpr2 = (1 - r - b) * (1 + b + r)
    sqarea = sqarea_triangle(one(T), r, b)
    k2 = max(zero(T), onembmr2 * fourbrinv)
    k = sqrt(k2)
    omemr2 = (1 - r) * (1 + r)
    onemr2mb2 = omemr2 - b2
    onemr2pb2 = omemr2 + b2
    if k2 ≥ 1
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
    kc = sqrt(abs(kc2))

    ## Compute uniform term

    (s0, kap0, kite_area2, kck), ∇s0 = compute_uniform_grad(b, r; b2, r2, sqarea, fourbrinv)

    flux = ld.g_n[begin] * s0
    ∇flux = ld.g_n[begin] * ∇s0
    dfdg[begin] = s0

    if ld.n_max == 0
        dfdg[begin] -= flux * ld.norm * π
        return flux * ld.norm, dfdg * ld.norm, ∇flux * ld.norm
    end

    ## compute linear term
    (s1, Eofk, Em1mKdm), ∇s1 = compute_linear_grad(b, r; k2, kc, kc2, r2, b2, br, fourbr, sqbr, onembmr2, onembpr2, onembmr2inv, sqonembmr2)

    flux += ld.g_n[begin + 1] * s1
    ∇flux = ∇flux + ld.g_n[begin + 1] * ∇s1
    dfdg[begin + 1] = s1

    if ld.n_max == 1
        dfdg[begin] -= flux * ld.norm * π
        dfdg[begin + 1] -= flux * ld.norm * π * 2/3
        return flux * ld.norm, dfdg * ld.norm, ∇flux * ld.norm
    end

    ## calculate quadratic term
    s2, ∇s2 = compute_quadratic_grad(b, r; s0, r2, b2, kap0, kite_area2, k2, ∇s0)
    flux += ld.g_n[begin + 2] * s2
    ∇flux = ∇flux + ld.g_n[begin + 2] * ∇s2
    dfdg[begin + 2] = s2
    if ld.n_max == 2
        dfdg[begin] -= flux * ld.norm * π
        dfdg[begin + 1] -= flux * ld.norm * π * 2/3
        return flux * ld.norm, dfdg * ld.norm, ∇flux * ld.norm
    end

    ## higher order (N > 3) terms require solving Mn and Nn integrals
    if k2 < 0.5 && ld.n_max > 3
        downwardM!(ld.Mn, ld.Mn_coeff; n_max=ld.n_max, sqonembmr2, sqbr, onemr2mb2, k, k2, sqarea, kap0, Eofk, Em1mKdm, kite_area2)
        if b < bcut
            Nn_lower!(ld.Nn, ld.Nn_coeff; Mn=ld.Mn, n_max=ld.n_max, kap0, kc, k2, sqbr, k, Eofk, Em1mKdm)
        end
    else
        upwardM!(ld.Mn; sqbr, n_max=ld.n_max, sqonembmr2, onemr2mb2, sqarea, k2, kap0, Eofk, Em1mKdm, kite_area2)
        if b < bcut
            Nn_raise!(ld.Nn, ld.Mn; n_max=ld.n_max, onembpr2, kap0, kc, k2, sqbr, k, Eofk, Em1mKdm)
        end
    end

    # compute remaining terms
    for n in 3:ld.n_max
        pofgn_M = 2 * r2 * ld.Mn[begin + n] - n / (n + 2) *
                       (onemr2mb2 * ld.Mn[begin + n] + sqarea * ld.Mn[begin + n - 2])
        flux -= ld.g_n[begin + n] * pofgn_M
        dfdg[begin + n] = -pofgn_M
        dpdr_M = 2 * r * ((n + 2) * ld.Mn[begin + n] - n * ld.Mn[begin + n - 2])

        if b < bcut
            dpdb_M = n * (ld.Mn[begin + n - 2] * (2 * r^3 + b^3 - b - 3 * r2 * b) + b * ld.Mn[begin + n] - 4 * r^3 * ld.Nn[begin + n - 2])
        else
            dpdb_M = n/b * ((ld.Mn[begin + n] - ld.Mn[begin + n - 2]) * (r2 + b2) + (b2 - r2)^2 * ld.Mn[begin + n - 2])
        end
        ∇flux = ∇flux - ld.g_n[begin + n] * SA[dpdb_M, dpdr_M]

    end
    dfdg[begin] -= flux * ld.norm * π
    dfdg[begin + 1] -= flux * ld.norm * π * 2/3
    return flux * ld.norm, dfdg * ld.norm, ∇flux * ld.norm
end

function Nn_lower!(Nn::AbstractVector, Nn_coeff; Mn, onembpr2, n_max, k, k2, kwargs...)
    Nn_series!(Nn, Nn_coeff; k2, sqonembmr2, n_max, k)

    for m in n_max-2:-1:2
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
        for n in 1:size(Nn_coeff, 3) - 1
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
    @show "raise"
    Nn_two!(Nn; kwargs...)
    for m in 2:n_max
        Nn[begin + m] = (Mn[begin + m] + m * onembpr2 * Nn[begin + m - 2]) / (m + 2)
    end
    return Nn
end

function Nn_two!(Nn::AbstractVector; kap0, kc, k2, sqbr, k, Eofk, Em1mKdm)
    if k2 ≤ 1
        Nn[begin] = 0.5 * kap0 - k * kc
        Nn[begin + 1] = 4/3 * sqbr * k2 * (2 * Em1mKdm - Eofk)
    else
        Nn[begin] = 0.5 * π
        Nn[begin + 1] = 4/3 * sqbr * k * (2 * Eofk - Em1mKdm)
    end

    return Nn
end