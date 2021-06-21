
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

    @test sprint(show, orbit) == "SimpleOrbit(P=3.0, T=1.0, t0=0.0, b=0.0)"

    @test sprint(show, "text/plain", orbit) == """
    SimpleOrbit
     period: 3.0
     duration: 1.0
     t0: 0.0
     b: 0.0"""
end

@testset "KeplerianOrbit" begin
    orbit = KeplerianOrbit(
        ρₛ = 2.0,
        Rₛ = 0.5,
        period = 2.0,
        ecc = 0.0,
        t₀ = 0.0,
        incl = π / 2.0,
        Ω = 0.0,
        ω = 0.0,
    )

    @test sprint(show, "text/plain", orbit) === """
    Keplerian Orbit
     a: 6.783710833739071 R⊙
     aRₛ: 13.567421667478142
     b: 8.307649758879776e-16
     ecc: 0.0
     P: 2.0 d
     ρₚ: 0.0 M⊙/R⊙³
     ρₛ: 11.810543837929684 g/cm³
     Rₚ: 0.001 R⊙
     Rₛ: 0.5 R⊙
     t₀: 0.0 d
     tₚ: -0.5 d
     t_ref: -0.5 d
     incl: 1.5707963267948966 rad
     Ω: 0.0 rad
     ω: 0.0 rad
     Mₛ: 1.0471975511965976 M⊙
     aₛ: 0.0 R⊙
     Mₚ: 0.0 M⊙
     aₚ: -6.783710833739071 R⊙"""
end
