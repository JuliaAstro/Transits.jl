using Transits: cel

# tests numerical accuracy
function test_cel!(kc, p, a, b)
    k2 = 1 - kc^2
    ell2 = collect(cel(k2, kc, p, a[1], a[2], a[3], b[1], b[2], b[3]))
    nϕ = 10^5
    if p < 0
        p0 = sqrt((kc^2 - p) / (1 - p))
        a0 = (a[1] - b[1]) / (1 - p)
        b0 = -(b[1] - a[1] * p) / (1 - p)^2 * (1 - kc^2) / p0 + a0 * p0
        a[1] = a0
        b[1] = b0 * p0
        p = p0^2
    end
    dϕ = 0.5 * π / nϕ
    ϕ = range(0.5 * dϕ, 0.5 * (π - dϕ); length=nϕ)

    ell1 = [
        cel(k2, kc, p, a[1], b[1]), cel(k2, kc, 1, a[2], b[2]), cel(k2, kc, 1, a[3], b[3])
    ]

    sϕ2 = @. sin(ϕ)^2
    cϕ2 = @. cos(ϕ)^2
    den = @. dϕ / sqrt(cϕ2 + kc^2 * sϕ2)
    ell3 = [
        sum((a[1] * cϕ2 + b[1] * sϕ2) ./ (cϕ2 + p * sϕ2) .* den),
        sum((a[2] * cϕ2 + b[2] * sϕ2) .* den),
        sum((a[3] * cϕ2 + b[3] * sϕ2) .* den),
    ]

    @test ell1 ≈ ell2 ≈ ell3
end

@testset "elliptic integrals - cel" begin
    kc = rand(rng)
    @testset "test negative p values" begin
        p = -rand(rng)
        a = rand(rng, 3)
        b = rand(rng, 3)
        test_cel!(kc, p, a, b)
    end
    @testset "test positive p values" begin
        p = rand(rng)
        a = rand(rng, 3)
        b = rand(rng, 3)
        test_cel!(kc, p, a, b)
    end
    @testset "test small kc values" begin
        p = rand(rng)
        a = rand(rng, 3)
        b = rand(rng, 3)
        kc = 1e-4
        test_cel!(kc, p, a, b)
    end
    @testset "test large kc values" begin
        p = rand(rng)
        a = rand(rng, 3)
        b = rand(rng, 3)
        kc = 1 - 1e-4
        test_cel!(kc, p, a, b)
    end
end
