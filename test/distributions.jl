using Distributions
using Bijectors
using Transits: Kipping13Transform
using HypothesisTests

@testset "Kipping13" begin
    N = 10000
    @test length(Kipping13()) == 2
    samps = rand(rng, Kipping13(), N)
    u1, u2 = eachrow(samps)

    @test all(@.(u1 + u2 < 1)) # flux is always positive
    @test all(@.(u1 + 2 * u2 > 0)) # flux is strictly monotonic increasing
    @test all(@.(u1 > 0)) # u1 always positive (flux is always positive)
    @test loglikelihood(Kipping13(), samps) == 0.0

    # make sure qs are Uniform
    q1 = @. (u1 + u2)^2
    q2 = @. 0.5 * u1 / (u1 + u2)
    for q in (q1, q2)
        t = ApproximateOneSampleKSTest(q, Uniform(0, 1))
        @test t.δ < 0.05
    end
end

@testset "Kipping13Transform" begin
    N = 1000
    @test bijector(Kipping13()) isa Kipping13Transform

    𝔅 = Kipping13Transform()
    xs = rand(Kipping13(), N)
    for x in eachcol(xs)
        y = 𝔅(x)
        xi = inverse(𝔅)(y)
        @test xi ≈ x
    end

    ys = 2 .* randn(2, N)
    @test all(eachcol(ys)) do y
        u1, u2 = inverse(𝔅)(y)

        u1 + u2 ≤ 1 && u1 + 2 * u2 ≥ 0 && u1 ≥ 0
    end
end
