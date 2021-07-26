import ChainRulesCore: frule, rrule
using LinearAlgebra


function compute_grad(ld::QuadLimbDark, b::S, r) where S
    T = float(S)
    bcut = 1e-3
    dfdg1 = zero(T)
    dfdg2 = zero(T)
    dfdg3 = zero(T)

    ## check for trivial cases
    if b ≥ 1 + r || iszero(r)
        # completely unobscured
        return one(T), SA[dfdg1, dfdg2, dfdg3], zero(T), zero(T)
    elseif r ≥ 1 + b
        # completely obscured
        return zero(T), SA[dfdg1, dfdg2, dfdg3], zero(T), zero(T)
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
            dfdr = ld.g_n[begin] * facd
        else
            flux = ld.g_n[begin] * onemr2 + 2 / 3 * ld.g_n[begin + 1] * sqrt1mr2^3
            dfdr = ld.g_n[begin] * facd + ld.g_n[begin + 1] * facd * sqrt1mr2
        end
        fac = 2 * r2 * onemr2
        if ld.n_max > 1
            flux -= ld.g_n[begin + 2] * fac
            dfdr += ld.g_n[begin + 2] * facd * 2 * (onemr2 - r2)
        end
        dfdr *= π * ld.norm
        dfdg = SA[(onemr2 - flux) * π * ld.norm,
                  2 / 3 * (sqrt1mr2^3 - flux) * π * ld.norm,
                  -fac]
        return flux * π * ld.norm, dfdg * ld.norm, zero(T), dfdr
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
    dfdg1 = s0

    if ld.n_max == 0
        dfdg1 -= flux * ld.norm * π
        dfdb, dfdr = ∇flux * ld.norm
        dfdg = SA[dfdg1, dfdg2, dfdg3]
        return flux * ld.norm, dfdg * ld.norm, dfdb, dfdr
    end

    ## compute linear term
    (s1, Eofk, Em1mKdm), ∇s1 = compute_linear_grad(b, r; k2, kc, kc2, r2, b2, br, fourbr, sqbr, onembmr2, onembpr2, onembmr2inv, sqonembmr2)

    flux += ld.g_n[begin + 1] * s1
    ∇flux += ld.g_n[begin + 1] * ∇s1
    dfdg2 = s1

    if ld.n_max == 1
        dfdg1 -= flux * ld.norm * π
        dfdg2 -= flux * ld.norm * π * 2 / 3
        dfdb, dfdr = ∇flux * ld.norm
        dfdg = SA[dfdg1, dfdg2, dfdg3]
        return flux * ld.norm, dfdg * ld.norm, dfdb, dfdr
    end

    ## calculate quadratic term
    s2, ∇s2 = compute_quadratic_grad(b, r; s0, r2, b2, kap0, kite_area2, k2, ∇s0)
    flux += ld.g_n[begin + 2] * s2
    ∇flux += ld.g_n[begin + 2] * ∇s2
    dfdg3 = s2

    dfdg1 -= flux * ld.norm * π
    dfdg2 -= flux * ld.norm * π * 2 / 3
    dfdb, dfdr = ∇flux * ld.norm
    dfdg = SA[dfdg1, dfdg2, dfdg3]
    return flux * ld.norm, dfdg * ld.norm, dfdb, dfdr
end

####

function frule((_, Δld, Δb, Δr), ::typeof(compute), ld::LD, b, r) where {LD <: QuadLimbDark}
    f, dfdg, dfdb, dfdr = compute_grad(ld, b, r)
    ∂g_n = dot(dfdg, Δld.g_n)
    return f, ∂g_n + dfdb * Δb + dfdr * Δr
end

function rrule(::typeof(compute), ld::LD, b, r) where {LD <: QuadLimbDark}
    f, dfdg, dfdb, dfdr = compute_grad(ld, b, r)
    function compute_pullback(Δf)
        ∂ld = Tangent{LD}(g_n=dfdg * Δf)
        ∂b = dfdb * Δf
        ∂r = dfdr * Δf
        return NoTangent(), ∂ld, ∂b, ∂r
    end
    return f, compute_pullback
end

function frule((_, Δu_n), ::Type{<:QuadLimbDark}, u_n::AbstractVector{T}) where T
    Ω = QuadLimbDark(u_n)
    ∇g_n = SA[-one(T) -1.5
               one(T)  2.0
              zero(T) -0.25]
    
    # calculate flux normalization factor, which only depends on first two terms
    N = length(u_n)
    if N == 0
        ∂g_n = @SVector zeros(T, 3)
    elseif N == 1
        ∂g_n = SA[-Δu_n[begin], Δu_n[begin], zero(T)]
    else
        Δu_n_full = SA[Δu_n[begin], Δu_n[begin + 1]]
        ∂g_n = ∇g_n * Δu_n_full
    end
    ∂Ω = Tangent{typeof(Ω)}(g_n=∂g_n)
    return Ω, ∂Ω
end


function rrule(::Type{<:QuadLimbDark}, u_n::AbstractVector{T}; maxiter=100) where T
    Ω = QuadLimbDark(u_n)
    ∇g_n = SA[-one(T)  one(T) zero(T)
                 -1.5     2.0  -0.25]
    N = length(u_n)
    function QuadLimbDark_pullback(Δld)
        ∂u = @view(∇g_n[begin:begin+N-1, :]) * Δld.g_n
        return NoTangent(), ∂u
    end
    return Ω, QuadLimbDark_pullback
end

