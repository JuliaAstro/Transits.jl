
u = [0.4, 0.26]

@testset "PolynomialLimbDark" begin
    ld = PolynomialLimbDark(u)

    @test repr("text/plain", ld) == """
    PolynomialLimbDark
     u_n: [-1.0, 0.4, 0.26]"""

    @test repr("text/html", ld) == """
    PolynomialLimbDark
     u<sup>n</sup>: [-1.0, 0.4, 0.26]"""

    @test sprint(show, ld) == "PolynomialLimbDark([-1.0, 0.4, 0.26])"

end

@testset "QuadLimbDark" begin
    ld2 = QuadLimbDark(u)
    @test repr("text/plain", ld2) == """
    QuadLimbDark
     u_n: [-1.0, 0.4, 0.26]"""

    @test repr("text/html", ld2) == """
    QuadLimbDark
     u<sup>n</sup>: [-1.0, 0.4, 0.26]"""

    @test sprint(show, ld2) == "QuadLimbDark([-1.0, 0.4, 0.26])"
end

@testset "SecondaryLimbDark" begin


end