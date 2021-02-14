using Transits.Orbits: relative_position, in_transit

@testset "simple orbit" begin
    period = 3.456
    t0 = 1.45
    b = 0.5
    duration = 0.12
    r_star = 1.345
    t = t0 .+ range(-2period, 2period, length=5000)
    m0 = @. abs((t - t0 + 0.5 * period) % period - 0.5 * period) < 0.5 * duration
    
    orbit = SimpleOrbit(;period, t0, b, duration, r_star)
    pos = @. relative_position(orbit, t)
    bs = map(v -> sqrt(v[1]^2 + v[2]^2), pos)
    zs = map(v->v[3], pos)
    m = @. vs ≤ r_star & zs > 0
    @test orbit.ref_time == -0.5
    @test m ≈ m0

    idxs = @. in_transit(orbit, t)
    @test all(bs[idxs] .≤ r_star)
    @test all(zs[idxs] .> 0)

    pos  = @. star_position(orbit, t)
    for vec in pos
        @test all(vec .≈ 0)
    end
end

@testset "simple limbdark" begin
    period = 3.456
    t0 = 1.45
    b = 0.5
    duration = 0.12
    r_star = 1.345


    t = t0 .+ range(-2period, 2period, length=5000)
    m0 = @. abs((t - t0 + 0.5 * period) % period - 0.5 * period) < 0.5 * duration
    orbit = SimpleOrbit(;period, t0, b, duration, r_star)


    ld = PolynomialLimbDark([0.2, 0.3])
    lc1 = ld.(orbit, t, 0.01)
    @test all(lc1[.!iszero.(m0)] ≈ 0.0)

    ldt = IntegratedLimbDark(ld)

    lc2 = ld.(orbit, t, 0.01, texp=0.01)
    @assert lc1 ≉ lc2
end