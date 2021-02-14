
function transit_poly_integrate(r, b, u, N)
    n = length(u)
    if r ≤ b - 1
        return 1.0
    end
    s_1 = max(0, b - r)
    s_2 = min(1, b + r)
    ds = (s_2 - s_1) / N
    s = range(s_1 + 0.5 * ds, s_2 - 0.5 * ds, length=N)
    
    # normalization
    norm = mapreduce(i -> 2 * u[i] / (2 + 3 * i + i^2), -, 1:n, init=1.0)

    flux = sum(s) do s_j
        if s_j < r - b
            dϕ = 2π
        else
            f = (s_j^2 + b^2 - r^2) / (2 * b * s_j)
            if abs(f) ≈ 1 && abs(f) > 1
                f = clamp(f, -1, 1)
            end
            dϕ = 2 * acos(f)
        end
        z = sqrt(1 - s_j^2)
        iμ = mapreduce(i -> u[i] * (1 - z)^i, -, 1:n, init=1.0)
        
        return s_j * ds * dϕ * iμ
    end

    flux /= π * norm
    flux = 1 - flux
    return flux
end

@testset "polynomial limb darkening" begin
    r0 = [
        fill(0.01, 5)
        fill(0.1, 5)
        fill(1, 3)
        fill(10, 3)
        fill(100, 3)
    ]
    b0 = [
        0, 0.1, 0.99, 1, 1.01,
        0, 0.1, 0.9, 1, 1.1,
        0, 1, 2,
        9, 10, 11,
        99, 100, 101
    ]
    n = 2 + ceil(Int, rand(rng) * 20)
    u = rand(rng, n)
    u *= rand(rng) / sum(u)

    N = 5000

    ld = PolynomialLimbDark(u)
    flux = @. ld(b0, r0)
    f_num = @. transit_poly_integrate(r0, b0, (u,), N)

    @info "Maximum difference of lightcurve for N=$N: $(maximum(abs,flux-f_num)))"
    @test flux ≈ f_num atol=1e-5

    # test vs. quad limb dark
    ld_ = PolynomialLimbDark(u[1:2])
    ld_q = QuadLimbDark(u[1:2])
    @test ld_.(b0, r0) == ld_q.(b0, r0)

    npts = 1000
    iu = 30
    u_n = ones(iu) / iu
    b0 = @. 3 * ((1:npts) / npts) - 0.5
    ld = PolynomialLimbDark(u_n)
    lc_ana = @. ld(abs(b0), 0.1)
    lc_num = @. transit_poly_integrate(0.1, abs(b0), (u_n,), N)

    @info "Maximum difference of lightcurve for N=$iu: $(maximum(abs,lc_ana-lc_num)))"
    @test lc_ana ≈ lc_num atol=1e-5

    

end


@testset "QuadLimbDark interface" begin
    # test truncation
    ld = (@test_logs (:warn, "Higher-order terms will be ignored") QuadLimbDark(ones(10)))
    @test ld.n_max == 2
    @test ld.g_n == QuadLimbDark(ones(2)).g_n
end

@testset "IntegratedLimbDark interface" for basis in [:legendre, :radau, :lobatto]
    u = [0.4, 0.26]
    ld = PolynomialLimbDark(u)
    ldt = IntegratedLimbDark(ld, basis=basis)

    @test IntegratedLimbDark(u).driver isa PolynomialLimbDark

    orbit = SimpleOrbit(period=3, duration=1)
    # sharp edges
    ts = [-0.5, -0.25, 0, 0.25, 0.5]
    above_or_below = [-1, 1, 1, 1, -1]

    # test no exposure time is exactly same
    @test ld.(orbit, ts, 0.01) == ldt.(orbit, ts, 0.01)

    # test exposure "smears" signal causing it to be
    # lower or higher then the exact signal
    f = @. ldt(orbit, ts, 0.01, [0.3 0.5])
    base = @. ld(orbit, ts, 0.01)
    @test sign.(f[:, 1] .- base) == above_or_below
    @test sign.(f[:, 2] .- base) == above_or_below
    # test higher exposure time is more smearing
    @test sign.(f[:, 2] .- f[:, 1]) == above_or_below

end

@testset "secondary" begin
    u = [0.4, 0.26]
    d1 = PolynomialLimbDark(u)
    d2 = PolynomialLimbDark(ones(2))
    ld = SecondaryLimbDark(d1, d2; brightness_ratio=0.1)
    ld2 = SecondaryLimbDark(u, ones(2))
    @test ld.primary_driver.g_n == ld2.primary_driver.g_n == d1.g_n
    @test ld.secondary_driver.g_n == ld2.secondary_driver.g_n == d2.g_n

    orbit = SimpleOrbit(period=4, duration=1)
    r = 0.01
    f1 = ld(orbit, 0, r)
    f2 = ld(orbit, 2, r)
    @test f1 ≈ f2 rtol=1e-3
end

@testset "reciprocity" begin
    u = [0.4, 0.26]
    d1 = PolynomialLimbDark(u)
    d2 = PolynomialLimbDark(ones(2))
    ld1 = IntegratedLimbDark(SecondaryLimbDark(d1, d2))
    ld2 = SecondaryLimbDark(IntegratedLimbDark(d1), IntegratedLimbDark(d2))

    orbit = SimpleOrbit(period=4, duration=1)
    r = 0.01
    t = range(-1.5, 1.5, length=1000)
    for texp in [nothing, 0.1, 0.3]
        @test ld1.(orbit, t, r; texp) ≈ ld2.(orbit, t, r; texp)
    end
end