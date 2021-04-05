using ForwardDiff
using FiniteDifferences
using Transits: compute_grad

# @testset "PolynomialLimbDark" begin
#     u = [0.4, 0.26, 0.1, 0.1]
#     test_frule(PolynomialLimbDark, u)
#     test_rrule(PolynomialLimbDark, u)

#     ld = PolynomialLimbDark(u)

#     test_frule(compute, ld, 0.0, 0.1)
#     test_rrule(compute, ld, 0.0, 0.1)
# end

# @testset "QuadLimbDark" begin
#     u = [0.4, 0.26]
#     test_frule(QuadLimbDark, u)
#     test_rrule(QuadLimbDark, u)

#     ld = QuadLimbDark(u)

#     test_frule(compute, ld, 0.0, 0.1)
#     test_rrule(compute, ld, 0.0, 0.1)
# end

function randbr(rng, condition)
    b = 2 * rand(rng)
    r = 2 * rand(rng)
    while !condition(b, r)
        b = 2 * rand(rng)
        r = 2 * rand(rng)
    end
    return b, r
end

function test_compute_grad(b, r, u_n)
    ld = PolynomialLimbDark(u_n)
    primal(X) = compute(ld, X[1], X[2])
    X = [b, r]
    grad_ForwardDiff = ForwardDiff.gradient(primal, X)
    _, _, dfdb, dfdr = compute_grad(ld, b, r)
    grad_analytical = [dfdb, dfdr]
    # grad_FiniteDiff = FiniteDifferences.grad(central_fdm(5, 1), primal, X)[1]

    # @test grad_analytical ≈ grad_FiniteDiff atol=1e-7
    @test grad_analytical ≈ grad_ForwardDiff atol=1e-8
end

@testset "`compute` grads" begin
    u_n = [0.4, 0.26]
    # b + r < 1
    b, r = randbr(rng, (b, r) -> b + r < 1)
    test_compute_grad(b, r, u_n)
    # b + r > 1
    b, r = randbr(rng, (b, r) -> b + r > 1)
    test_compute_grad(b, r, u_n)

    # special cases
    test_compute_grad(0.5, 0.5, u_n) # r = b = 1/2
    test_compute_grad(0.3, 0.0, u_n) # b = 0
    test_compute_grad(0.2, 0.8, u_n) # r + b = 1
    test_compute_grad(0.3, 0.3, u_n) # r = b < 1/2
    test_compute_grad(3.0, 3.0, u_n) # r = b > 1/2 
end