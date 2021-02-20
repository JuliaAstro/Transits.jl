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

function compute_uniform_grad(b::T, r; kap0, kite_area2, kwargs...) where T
    s0 = compute_uniform(b, r; kap0, kite_area2, kwargs...)
    if b ≤ 1 - r
        ∇s0 = SA[-2 * π * r, zero(T)]
    else
        ∇s0 = SA[-2 * r * kap0, kite_area2 / b]
    end
    return s0, ∇s0
end

function frule((_, _, _), ::typeof(compute_uniform), b, r; kwargs...)
    compute_uniform_grad(b, r; kwargs...)
end


function compute_linear_grad(b::T, r; k2, kc, kc2, onembmr2, onembpr2, onembmr2inv, r2=r^2, b2=b^2, br=b*r, fourbr=4*br, sqbr=sqrt(br)) where T
    if iszero(b) # case 10
        Λ1 = -2 * π * sqrt(1 - r2)^3
        Eofk = 0.5 * π
        Em1mKdm = 0.25 * π

        ∇s1 = SA[-2 * π * r * sqrt1mr2, zero(T)]
    elseif b == r
        if r == 0.5 # case 6
            Λ1 = π - 4 / 3 - 2 * (b - 0.5) + 6 * (r - 0.5)
            Eofk = one(T)
            Em1mKdm = one(T)
            ∇s1 = SA[-2, 2/3]
        elseif r < 0.5 # case 5
            m = 4 * r2
            Eofk = cel(m, one(T), one(T), 1 - m)
            Em1mKdm = cel(m, one(T), one(T), zero(T))
            Λ1 = π + 2 / 3 * ((2 * m - 3) * Eofk - m * Em1mKdm) +
                 (b - r) * 4 * r * (Eofk - 2 * Em1mKdm)
            ∇s1 = SA[-4 * r * Eofk, -4/3 * r * (Eofk -2 * Em1mKdm)]
        else # case 7
            m = 4 * r2
            minv = inv(m)
            Eofk = cel(minv, one(T), one(T), 1 - minv)
            Em1mKdm = cel(minv, one(T), one(T), zero(T))
            Λ1 = π + 1 / 3 * ((2 * m - 3) * Em1mKdm - m * Eofk) / r -
                 (b - r) * 2 * (2 * Eofk - Em1mKdm)
            ∇s1 = SA[-2 * Em1mKdm, 2/3 * (2 * Eofk - Em1mKdm)]
        end
    else
        if b + r > 1 # case 2, case 8
            Πofk, Eofk, Em1mKdm = cel(k2, kc, (b - r)^2 * kc2, zero(T), one(T), one(T), 3 * kc2 * (b - r) * (b + r), kc2, zero(T))
            Λ1 = onembmr2 * (Πofk + (-3 + 6 * r2 + 2 * br) * Em1mKdm - fourbr * Eofk) / (3 * sqbr)
            ∇s1 = SA[-2 * r * onembmr2 * Em1mKdm * sqbrinv, 2 * r * onembmr2 * (-Em1mKdm + 2 * Eofk) * sqbrinv/3]
        elseif b + r < 1 # case 3 case 9
            bmrdbpr = (b - r) / (b + r)
            μ = 3 * bmrdbpr * onembmr2inv
            p = bmrdbpr^2 * onembpr2 * onembmr2inv
            Πofk, Eofk, Em1mKdm = cel(inv(k2), kc, p, 1 + μ, one(T), one(T), p + μ, kc2, zero(T))
            Λ1 = 2 * sqrt(onembmr2) * (onembpr2 * Πofk - (4 - 7 * r2 - b2) * Eofk) / 3
            ∇s1 = SA[-4 * r * sqonembmr2 * Eofk, -4/3 * r * sqonembmr2 * (Eofk - 2 * Em1mKdm)]
        else # case
            Λ1 = 2 * acos(1 - 2 * r) - 2 * π * (r > 0.5) - (4/3 * (3 + 2 * r - 8 * r2) + 8 * (r + b - 1) * r) * sqrt(r * (1 - r))
            Eofk = one(T)
            Em1mKdm = one(T)
            g = -8 * r * sqrt(r * (1 - r))
            ∇s1 = SA[g, -g/3]
        end
    end

    flux = ((1 - T(r > b)) * 2π - Λ1) / 3

    return (flux, Eofk, Em1mKdm), ∇s1
end

function frule((_, _, _), ::typeof(compute_linear), b, r; kwargs...)
    compute_linear_grad(b, r; kwargs...)
end


function compute_quadratic_grad(b, r; s0, r2, b2, kap0, kite_area2, k2, ∇s1)
    r2pb2 = r2 + b2
    η2 = r2 * (r2pb2 + b2)
    if k2 > 1
        four_pi_eta = 2 * π * (η2 - 1)
        ∇η = 4 * π * SA[2 * r * r2pb2, 2 * b * r2]
    else
        Πmkap1 = atan(kite_area2, (r - 1) * (r + 1) - b2)
        four_pi_eta = 2 * (-Πmkap1 + η2 * kap0 - 0.25 * kite_area2 * (1 + 5 * r2 + b2))
        ∇η = SA[8 * r * (r2pb2 * kap0 - kite_area2), 
                2/b * (4 * b2 * r2 * kap0 - (1 + r2pb2) * kite_area2)]
    end
    s2 = 2 * s0 + four_pi_eta
    ∇s2 = 2 * ∇s1 + ∇η
    return  s2, ∇s2
end

function frule((_, _, _), ::typeof(compute_quadratic), b, r; kwargs...)
    compute_quadratic_grad(b, r; kwargs...)
end
