import ChainRulesCore: frule, rrule


function frule((_, _), ::typeof(compute_gn), u_n::AbstractVector{T}) where T
    compute_gn_jacobian(u_n)
end

function compute_gn_jacobian(u_n::AbstractVector{T}) where T
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

function compute_gn_jacobian(u_n::StaticVector{3,T}) where T
    g_n = compute_gn(u_n)
    Δg_n = SA[zero(T) -one(T)  -1.5
              zero(T)  one(T)   2.0 ; 
              zero(T)  zero(T) -0.25]
    return g_n, Δg_n
end
