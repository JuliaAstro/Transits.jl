
@testset "Kipping13" begin
    N = 10000
    @test length(Kipping13()) = 2
    u1, u2 = rand(rng, Kipping13(), N) |> eachrow

    @test all(u1 .+ u2 .< 1) # flux is always positive
    @test all(u1 .+ 2 .* u2 .> 0) # flux is strictly monotonic increasing
    @test all(u1 .> 0) # u1 always positive (flux is always positive)
end
