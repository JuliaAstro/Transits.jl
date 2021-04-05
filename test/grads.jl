using ForwardDiff
using FiniteDifferences
using Transits: compute_grad

# @testset "PolynomialLimbDark" begin
#     u = [0.4, 0.26, 0.1, 0.1]
#     # test_frule(PolynomialLimbDark, u)
#     # test_rrule(PolynomialLimbDark, u)

#     ld = PolynomialLimbDark(u)

#     test_frule(compute, ld, 0.0, 0.1; fdm=forward_fdm(10, 1))
#     test_rrule(compute, ld, 0.0, 0.1; fdm=forward_fdm(10, 1))
# end

# @testset "QuadLimbDark" begin
#     u = [0.4, 0.26]
#     test_frule(QuadLimbDark, u)
#     test_rrule(QuadLimbDark, u)

#     ld = QuadLimbDark(u)

#     test_frule(compute, ld, 0.0, 0.1)
#     test_rrule(compute, ld, 0.0, 0.1)
# end

function finite_diff(b_, r_, u_n_, law=PolynomialLimbDark; diff=big(1e-18))
    b = big(b_)
    r = big(r_)
    u_n = big.(u_n_)
    ld = law(u_n)
    f = ld(b, r)

    # modulate r
    f₊ = ld(b, r + diff)
    if r == 1 && length(u_n) > 2
        dfdr = (f₊ - f) / diff
    else
        f₋ = ld(b, r - diff)
        dfdr = (f₊ - f₋) / (2 * diff)
    end

    # modulate b
    f₊ = ld(b + diff, r)
    if b > diff
        f₋ = ld(b - diff, r)
        dfdb = (f₊ - f₋) / (2 * diff)
    elseif iszero(b)
        dfdb = zero(b)
    else
        dfdb = (f₊ - f) / diff
    end
    
    return Float64.([dfdb, dfdr])
end

function test_compute_grad(b, r, u_n, law=PolynomialLimbDark)
    ld = law(u_n)
    _, _, dfdb, dfdr = @inferred compute_grad(ld, b, r)
    grad_analytical = [dfdb, dfdr]

    atol = r != 1 ? 1e-7 : 3e-6

    @testset "compute_grad - b=$b, r=$r" begin
        # test against BigFloat for numerical accuracy
        _, _, dfdb_big, dfdr_big = compute_grad(ld, big(b), big(r))
        grad_big = [dfdb_big, dfdr_big]
        @test grad_analytical ≈ grad_big atol = atol

        # test against finite diff for correctness
        grad_numerical = finite_diff(b, r, u_n, law)
        @test grad_analytical ≈ grad_numerical atol = atol
    end
end

logspace(start, stop, N) = 10 .^range(log10(start), log10(stop), length=N)

test_coeffs = [[0.0], [1.0], [2.0, -1.0], [3.0, -3.0, 1.0], fill(0.1, 4), fill(0.1, 5), fill(0.1, 10)]
test_names = ["uniform", "linear", "quadratic", "cubic", "quartic", "quintic", "10th order"]

@testset "`compute` grads - $name" for (name, u_n) in zip(test_names, test_coeffs)
    # b + r < 1
    rs = [0.01, 0.1, 0.5, 0.999999]
    nb = 50
    ϵ, δ = 1e-9, 1e-3
    for r in rs
        bs = abs.([
            range(0, ϵ, length=nb);
            range(ϵ, δ, length=nb);
            range(δ, r - δ, length=nb);
            -logspace(δ, ϵ, nb) .+ r;
            range(r - ϵ, r + ϵ, length=nb);
            logspace(ϵ, δ, nb) .+ r;
            range(r + δ, 1 - r - δ, length=nb);
            -logspace(δ, ϵ, nb) .+ (1 - r);
            range(1 - r - ϵ, 1 - r + ϵ, length=nb);
            logspace(ϵ, δ, nb) .+ (1 - r);
            range(1 - r + δ, 1 + r - δ, length=nb);
            -logspace(δ, ϵ, nb) .+ (1 + r);
            range(1 + r - ϵ, 1 + r - 1e-13, length=nb)
        ])
        test_compute_grad.(bs, r, (u_n,), PolynomialLimbDark)
        if length(u_n) < 3
            test_compute_grad.(bs, r, (u_n,), QuadLimbDark)
        end
    end

    rs = [1, 1.000001, 2, 10, 100]
    for r in rs
        bs = abs.([
            r - 1 + 1e-13;
            logspace(ϵ, δ, nb) .+ (r - 1);
            range(r - 1 + δ, r - δ, length=nb);
            -logspace(δ, ϵ, nb) .+ r;
            range(r - ϵ, r + ϵ, length=nb);
            logspace(ϵ, δ, nb) .+ r;
            range(r + δ, r + 1 - δ, length=nb);
            -logspace(δ, ϵ, nb) .+ (r + 1);
            r + 1 - 1e-13
        ])
        test_compute_grad.(bs, r, (u_n,), PolynomialLimbDark)
        if length(u_n) < 3
            test_compute_grad.(bs, r, (u_n,), QuadLimbDark)
        end
    end
end
