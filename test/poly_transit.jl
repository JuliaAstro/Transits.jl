
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
