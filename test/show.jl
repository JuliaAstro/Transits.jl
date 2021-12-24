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

    @test sprint(show, orbit) == "SimpleOrbit(P=3, T=1, t0=0, b=0)"

    @test sprint(show, "text/plain", orbit) == """
    SimpleOrbit
     period: 3
     duration: 1
     t0: 0
     b: 0"""
end

@testset "KeplerianOrbit" begin
    orbit = KeplerianOrbit(
        ρ_star = 1.0,
        R_star = cbrt(3.0/(4.0*π)),
        period = 2.0,
        ecc = 0.0,
        t_0 = 0.0,
        incl = π / 2.0,
        Ω = 0.0,
        ω = 0.0,
    )

    orbit_unit = KeplerianOrbit(
        ρ_star = 1.0u"g/cm^3",
        R_star = cbrt(3.0 / (4.0 * π))u"Rsun",
        period = 2.0u"d",
        ecc = 0.0,
        t_0 = 0.0u"d",
        incl = π / 2.0,
        Ω = 0.0,
        ω = 0.0,
    )

    @test repr("text/plain", orbit) === """
    Keplerian Orbit
     P:      2.0000 d
     t₀:     0.0000 d
     tₚ:     -0.5000 d
     t_ref:  -0.5000 d
     τ:      nothing
     a:      6.6802 R⊙
     aₚ:     -6.6802 R⊙
     aₛ:     0.0000 R⊙
     Rₚ:     0.0000 R⊙
     Rₛ:     0.6204 R⊙
     ρₚ:     NaN M⊙/R⊙³
     ρₛ:     1.0000 M⊙/R⊙³
     r:      0.0000
     aRₛ:    10.7685
     b:      0.0000
     ecc:    0.0000
     cos(i): 0.0000
     sin(i): 1.0000
     cos(ω): 1.0000
     sin(ω): 0.0000
     cos(Ω): 1.0000
     sin(Ω): 0.0000
     i:      1.5708 rad
     ω:      0.0000 rad
     Ω:      0.0000 rad
     Mₚ:     0.0000 M⊙
     Mₛ:     1.0000 M⊙"""

    @test repr("text/plain", orbit_unit) === """
    Keplerian Orbit
     P:      2.0000 d
     t₀:     0.0000 d
     tₚ:     -0.5000 d
     t_ref:  -0.5000 d
     τ:      nothing
     a:      3.6958 R⊙
     aₚ:     -3.6958 R⊙
     aₛ:     0.0000 R⊙
     Rₚ:     0.0000 R⊙
     Rₛ:     0.6204 R⊙
     ρₚ:     NaN M⊙ R⊙^-3
     ρₛ:     0.1693 M⊙ R⊙^-3
     r:      0.0000
     aRₛ:    5.9576
     b:      0.0000
     ecc:    0.0000
     cos(i): 0.0000
     sin(i): 1.0000
     cos(ω): 1.0000
     sin(ω): 0.0000
     cos(Ω): 1.0000
     sin(Ω): 0.0000
     i:      1.5708 rad
     ω:      0.0000 rad
     Ω:      0.0000 rad
     Mₚ:     0.0000 M⊙
     Mₛ:     0.1693 M⊙"""
end
