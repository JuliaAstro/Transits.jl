
u = [0.4, 0.26]

@testset "PolynomialLimbDark" begin
    ld = PolynomialLimbDark(u)

    @test sprint(show, ld) == "PolynomialLimbDark([-1.0, 0.4, 0.26])"

    @test repr("text/plain", ld) == """
    PolynomialLimbDark
     u_n: [-1.0, 0.4, 0.26]"""

    @test repr("text/html", ld) == """
    PolynomialLimbDark
     u<sub>n</sub>: [-1.0, 0.4, 0.26]"""
end

@testset "QuadLimbDark" begin
    ld = QuadLimbDark(u)

    @test sprint(show, ld) == "QuadLimbDark([-1.0, 0.4, 0.26])"

    @test sprint(show, "text/plain", ld) == """
    QuadLimbDark
     u_n: [-1.0, 0.4, 0.26]"""

    @test sprint(show, "text/html", ld) == """
    QuadLimbDark
     u<sub>n</sub>: [-1.0, 0.4, 0.26]"""
end

@testset "SecondaryLimbDark" begin
    ld = SecondaryLimbDark(u)

    @test sprint(show, ld) == "SecondaryLimbDark(PolynomialLimbDark, PolynomialLimbDark, 1.0)"

    @test sprint(show, "text/plain", ld) == """
    SecondaryLimbDark
     primary: PolynomialLimbDark([-1.0, 0.4, 0.26])
     secondary: PolynomialLimbDark([-1.0, 0.4, 0.26])
     ratio: 1.0"""
end

@testset "IntegratedLimbDark" begin
    ld = IntegratedLimbDark(u)

    @test sprint(show, ld) == "IntegratedLimbDark(PolynomialLimbDark, 21)"

    @test sprint(show, "text/plain", ld) == """
    IntegratedLimbDark
     driver: PolynomialLimbDark([-1.0, 0.4, 0.26])
     N: 21"""
end

@testset "SimpleOrbit" begin
    orbit = SimpleOrbit(duration=1, period=3)

    @test sprint(show, orbit) == "SimpleOrbit(P=3.0, T=1.0, t0=0.0, b=0)"

    @test sprint(show, "text/plain", orbit) == """
    SimpleOrbit
     period: 3.0
     duration: 1.0
     t0: 0.0
     b: 0"""
end
