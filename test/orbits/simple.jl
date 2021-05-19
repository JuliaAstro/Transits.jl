using Transits.Orbits: relative_position

@testset "simple orbit" begin
    period = 3.456
    t0 = 1.45
    b = 0.5
    duration = 0.12
    t = t0 .+ range(-2.0*period, 2.0*period, length=5000)
    m0 = @. abs(mod(t - t0 + 0.5*period, period) - 0.5*period) < 0.5 * duration

    orbit = SimpleOrbit(;period, t0, b, duration)
    pos = @. relative_position(orbit, t)
    bs = map(v -> sqrt(v[1]^2.0 + v[2]^2.0), pos)
    zs = map(v->v[3], pos)
    m = @. (bs ≤ 1.0) & (zs > 0.0)
    @test orbit.ref_time == -0.278
    @test m ≈ m0

    #idxs = in_transit(orbit, t)
    #@test all(bs[idxs] .≤ 1.0)
    #@test all(zs[idxs] .> 0.0)

    # TODO: Do we want this?
    #pos  = @. star_position(orbit, t)
    #for vec in pos
    #    @test all(vec .≈ 0.0)
    #end
end

@testset "simple limbdark" begin
    period = 3.456
    t0 = 1.45
    b = 0.5
    duration = 0.12

    t = t0 .+ range(-2.0*period, 2.0*period, length=5000)
    m0 = @. abs(mod(t - t0 + 0.5*period, period) - 0.5*period) < 0.5 * duration
    orbit = SimpleOrbit(;period, t0, b, duration)
    ld = PolynomialLimbDark([0.2, 0.3])

    lc1 = @. ld.(orbit, t, 0.01)
    @test all(lc1[.!m0] .== 1.0)

    ldt = IntegratedLimbDark(ld)
    lc2 = @. ldt.(orbit, t, 0.01, 0.01)
    @assert lc1 ≉ lc2
end
