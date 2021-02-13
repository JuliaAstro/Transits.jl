
function Mn_num(k2::T, m::Int64) where T
    f(x) = sqrt(k2 - sin(x)^2)^m
    if k2 < 1
        kap2 = T(asin(sqrt(big(k2))))
    else
        kap2 = 0.5 * π
    end
    Mn, err = quadgk(f, -kap2, kap2, rtol=1e-15)
    return Mn
end

function test_Mn(r::T, b::T) where T
    n_max = 30
    u = ones(n_max + 2)
    prec_frac = zeros(n_max + 3)
    prec_abs = zeros(n_max + 3)

    ld = PolynomialLimbDark(u)
    ld_b = PolynomialLimbDark(big.(u))

    t = ld(b, r)
    t_b = ld_b(b, r)


    onembpr2 = (1 - r - b) * (1 + b + r)
    k2 = onembpr2 * inv(4 * b * r) + 1

    k2 < 0 && return prec_frac, prec_abs
    
    sqbr = sqrt(b * r)

    for m in 0:ld.n_max
        Mnn = Mn_num(k2, m) * (2 * sqbr)^m
        Mn = ld.Mn[m + 1]
        Mn_b = Float64(ld_b.Mn[m + 1])
        t1 = isapprox(Mnn, Mn; atol=1e-15, rtol=1e-6)
        t2 = isapprox(Mn, Mn_b; atol=1e-15, rtol=1e-6)
        if !t1
            diff = Mnn - Mn
            @warn "Mn discrepancy" b r k2 Mnn Mn diff
        end
        if !t2
            diff = Mn - Mn_b
            @warn "Mn discrepancy" b r k2 Mn Mn_b diff
        end
        # @test Mnn ≈ Mn atol=1e-15 rtol=1e-6
        # @test Mn ≈ Mn_b atol=1e-15 rtol=1e-6
        @test t1 && t2
        prec_frac[m + 1] = Mn / Mn_b - 1
        prec_abs[m + 1] = asinh(Mn) - asinh(Mn_b)
    end

    return prec_frac, prec_abs
end


@testset "Mn compute" begin
    r = 100.0
    ϵ = 1e-8 
    bs = [r - 1 + ϵ, r, r + 1 - ϵ]
       
    prec_frac, prec_abs = test_Mn.(r, bs)
       

    r = 0.01
    bs = [ϵ, r, 1 - r - ϵ, 1 - r + ϵ, 1, r + 1 - ϵ]
    prec_frac, prec_abs = test_Mn.(r, bs)

    # test random points
    ntest = 10
    for i in 1:ntest
        r = 2.0
        b = 0.25
        while r > 1 + b || r < b - 1
            r = 2 * rand(rng)
            b = 2 * rand(rng)
        end
        k2 = (1-b+r)*(1+b-r)/(4*b*r)
        prec_frac, prec_abs = test_Mn(r, b)
    end
end
