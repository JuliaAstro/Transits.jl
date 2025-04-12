var documenterSearchIndex = {"docs":
[{"location":"gettingstarted/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"gettingstarted/#Usage","page":"Getting Started","title":"Usage","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"using Transits\n\norbit = SimpleOrbit(period=3, duration=1)\nu = [0.4, 0.26] # quad limb dark\nld = PolynomialLimbDark(u)\n\nt = range(-1, 1, length=1000) # days from t0\nrs = range(0, 0.2, length=10) # radius ratio\n\nfluxes = @. ld(orbit, t, rs')","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"gettingstarted/#Integrated-and-Secondary-Curves","page":"Getting Started","title":"Integrated and Secondary Curves","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"IntegratedLimbDark can be used to numerically integrate each light curve exposure in time","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"ld = IntegratedLimbDark([0.4, 0.26])\norbit = SimpleOrbit(period=3, duration=1)\nt = range(-1, 1, length=1000)\ntexp = [0.1 0.2 0.3]\n# no extra calculations made\nflux = @. ld(orbit, t, 0.2)\n# use quadrature to find time-averaged flux for each t\nflux_int = @. ld(orbit, t, 0.2, texp) ","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"SecondaryLimbDark can be used to generate secondary eclipses given a surface brightness ratio","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"ld = SecondaryLimbDark([0.4, 0.26], brightness_ratio=0.1)\nld_int = IntegratedLimbDark(ld) # composition works flawlessly\n\norbit = SimpleOrbit(period=4, duration=1)\nt = range(-1.25, 2.75, length=1000)\nrs = range(0.01, 0.1, length=6)\n\nf = @. ld(orbit, t, rs')\nf_int = @. ld_int(orbit, t, rs', texp=0.3)","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"gettingstarted/#Using-Units","page":"Getting Started","title":"Using Units","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"Units from Unitful.jl are a drop-in substitution for numbers","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"using Unitful\norbit = SimpleOrbit(period=10u\"d\", duration=5u\"hr\")\nt = range(-6, 6, length=1000)u\"hr\"\nflux = @. ld(orbit, t, 0.1)","category":"page"},{"location":"gettingstarted/#Gradients","page":"Getting Started","title":"Gradients","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"Gradients are provided in the form of chain rules. The easiest way to access them is using an automatic differentiation (AD) library like ForwardDiff.jl or Zygote.jl.","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"using Zygote\n\nts = range(-1, 1, length=1000) # days from t0\nror = 0.1\nu_n = [0.4, 0.26]\n\norbit = SimpleOrbit(period=3, duration=1)\nlightcurve(X) = compute(PolynomialLimbDark(X[3:end]), orbit, X[1], X[2])\n\n# use Zygote for gradient\nflux = [lightcurve([t, ror, u_n...]) for t in ts]\ngrads = mapreduce(hcat, ts) do t\n    grad = lightcurve'([t, ror, u_n...])\n    return grad === nothing ? zeros(4) : grad\nend","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"api/#API/Reference","page":"API/Reference","title":"API/Reference","text":"","category":"section"},{"location":"api/#Index","page":"API/Reference","title":"Index","text":"","category":"section"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"","category":"page"},{"location":"api/#Light-Curves","page":"API/Reference","title":"Light Curves","text":"","category":"section"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"AbstractLimbDark\n(AbstractLimbDark)(args...; kwargs...)\nPolynomialLimbDark\nQuadLimbDark\nIntegratedLimbDark\nSecondaryLimbDark\ncompute\ncompute(::AbstractLimbDark, ::Orbits.AbstractOrbit, t, r)","category":"page"},{"location":"api/#Transits.AbstractLimbDark","page":"API/Reference","title":"Transits.AbstractLimbDark","text":"AbstractLimbDark\n\nA limb dark law need only need to implement compute(::Law, b, r) to extend the limb darkening interface.\n\nSee also\n\ncompute\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.AbstractLimbDark-Tuple","page":"API/Reference","title":"Transits.AbstractLimbDark","text":"(::AbstractLimbDark)(b, r)\n\nAn alias for calling compute\n\nExamples\n\njulia> ld = PolynomialLimbDark([0.4, 0.26]);\n\njulia> ld(0, 0.01)\n0.9998785437247428\n\n\n\n\n\n","category":"method"},{"location":"api/#Transits.PolynomialLimbDark","page":"API/Reference","title":"Transits.PolynomialLimbDark","text":"PolynomialLimbDark(u::AbstractVector)\n\nPolynomial limb darkening using analytical integrals. The length of the u vector is equivalent to the order of polynomial used; e.g., [0.2, 0.3] corresponds to quadratic limb darkening.\n\nMathematical form\n\nI(mu) propto 1 - u_1(1-mu) - u_2(1-mu)^2 - dots - u_N(1-mu)^N\n\nwhich is equivalent to the series\n\nI(mu) propto -sum_i=0^Nu_i(1-mu)^i\n\nwith the definition u_0 equiv -1.\n\nExamples\n\nu = [0.4, 0.26] # quadratic and below is 100% analytical\nld = PolynomialLimbDark(u)\nld(0.1, 0.01)\n\n# output\n0.9998787880717668\n\nu2 = vcat(u, ones(12) ./ 12)\nld2 = PolynomialLimbDark(u2)\nld2(0.1, 0.01)\n\n# output\n0.9998740059086433\n\nReferences\n\nAgol, Luger, Foreman-Mackey (2020)\"Analytic Planetary Transit Light Curves and Derivatives for Stars with Polynomial Limb Darkening\"\n\nLuger et al. (2019)\"starry: Analytic Occultation Light Curves\"\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.QuadLimbDark","page":"API/Reference","title":"Transits.QuadLimbDark","text":"QuadLimbDark(u::AbstractVector)\n\nA specialized implementation of PolynomialLimbDark with a maximum of two terms (quadratic form). This has a completely closed-form solution without any numerical integration. This means there are no intermediate allocations and reduced numerical error.\n\nMathematical form\n\nI(mu) propto 1 - u_1(1-mu) - u_2(1-mu)^2\n\nwarning: Higher-order terms\nHigher-order terms will be ignored; no error will be thrown\n\nExamples\n\nld = QuadLimbDark(Float64[]) # constant term only\n\nb = [0, 1, 2] # impact parameter\nr = 0.01 # radius ratio\nld.(b, r)\n\n# output\n3-element Vector{Float64}:\n 0.9999\n 0.9999501061035608\n 1.0\n\nld = QuadLimbDark([0.4, 0.26]) # max two terms\nld.(b, r)\n\n# output\n3-element Vector{Float64}:\n 0.9998785437247428\n 0.999974726693709\n 1.0\n\nReferences\n\nSee references for PolynomialLimbDark\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.IntegratedLimbDark","page":"API/Reference","title":"Transits.IntegratedLimbDark","text":"IntegratedLimbDark(limbdark; N=21, basis=:legendre)\nIntegratedLimbDark(u; kwargs...)\n\nComputes the time-averaged flux in the middle of an exposure by wrapping a limb darkening law limbdark with a quadrature scheme. For each time step t, N extra points are super-sampled from t-texp/2 to t+texp/2and the time-averaged flux is calculated via quadrature.\n\nIf a set of limb darkening coefficients, u, is provided, a PolynomialLimbDark law will be used by default.\n\nMathematical form\n\nbarF(t) = frac1Delta tint_t-Delta t  2^t+Delta t  2F(t)dt\n\nwhere F is the wrapped limb darkening law and Delta t is the exposure time.\n\nQuadrature\n\nThe integration is approximated via Guassian quadrature\n\nfrac1Delta t intF(t)dt approx frac12sum_i^Nw_i * F(fracDelta t2xi_i + t)\n\nwhere the weights w_i and nodes ξ_i are defined by the given quadrature rule. The nodes are defined by evaluating orthogonal polynomials N times between -1 and 1. Notice the change of interval required to go from the natural bounds of the orthogonal polynomial basis, -1, 1, to the range defined by the exposure time.\n\nThe following bases are available from FastGaussQuadrature.jl. In addition, a function can be passed which calculates nodes, weights = f(N).\n\n:legendre - Legendre polynomial base on the open (-1, 1)\n:radau - Legendre polynomial base on the semi-open [-1, 1) interval\n:lobatto - Legendre polynomial base on the closed [-1, 1] interval\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.SecondaryLimbDark","page":"API/Reference","title":"Transits.SecondaryLimbDark","text":"SecondaryLimbDark(primary::AbstractLimbDark,\n                  secondary::AbstractLimbDark; \n                  brightness_ratio=1)\nSecondaryLimbDark(u_p::AbstractVector, u_s=u_p; kwargs...)\n\nCompose two limb darkening laws together to add a secondary eclipse. If vectors of coefficients are provided, laws will automatically be constructed using PolynomialLimbDark. The surface brightness ratio is given in terms of the host; e.g., if the companion is half as bright as the host, the ratio would be 0.5.\n\nnote: Interface\nSecondaryLimbDark only works with an orbit, since the companion's reference frame needs to be calculated. This means you can't call it using an impact parameter like ld(b, r) directly.\n\nMathematical form\n\nf(t r) = frac2f_p(t r) + eta r^2 f_s(t r)1 + f_p(t r) + eta r^2 f_s(t r)\n\nwhere f_p is to the primary flux, f_s is to the secondary flux, and eta is the surface brightness ratio. t and r correspond to the time and radius ratio from the companion's reference frame.\n\nExamples\n\n# equal size and limb darkening\nr = 1.0\nu = [0.4, 0.26]\n# companion is 1/10 as bright\nbrightness_ratio = 0.1\nld = SecondaryLimbDark(u; brightness_ratio)\norbit = SimpleOrbit(period=2, duration=0.5)\nfp = ld(orbit, 0, r) # primary egress\nfs = ld(orbit, 1, r) # secondary egress\n\nfp ≈ brightness_ratio * fs\n\n# output\ntrue\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.compute","page":"API/Reference","title":"Transits.compute","text":"compute(::AbstractLimbDark, b, r; kwargs...)\n\nCompute the relative flux for the given impact parameter b and radius ratio r. The impact parameter is unitless. The radius ratio is given in terms of the host; e.g., if the companion is half the size of the host, r=0.5.\n\n\n\n\n\n","category":"function"},{"location":"api/#Transits.compute-Tuple{AbstractLimbDark, Transits.Orbits.AbstractOrbit, Any, Any}","page":"API/Reference","title":"Transits.compute","text":"compute(::AbstractLimbDark, orbit::AbstractOrbit, t, r)\n\nCompute the relative flux by calculating the impact parameter at time t from the given orbit. The time needs to be compatible with the period of the orbit, nominally in days.\n\nExamples\n\njulia> ld = PolynomialLimbDark([0.4, 0.26]);\n\njulia> orbit = SimpleOrbit(period=3, duration=1);\n\njulia> ld(orbit, 0, 0.1) # primary egress\n0.9878664434953113\n\njulia> ld(orbit, 0.1, 0.1) # 0.1 d\n0.9879670695533511\n\nthis works effortlessly with libraries like Unitful.jl\n\njulia> using Unitful\n\njulia> orbit = SimpleOrbit(period=3u\"d\", duration=3u\"hr\");\n\njulia> ld(orbit, 0u\"d\", 0.1)\n0.9878664434953113\n\n\n\n\n\n","category":"method"},{"location":"api/#Gradients","page":"API/Reference","title":"Gradients","text":"","category":"section"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"Gradients and jacobians are integrated directly into ChainRules.jl via frules and rrules. For most users, this just means using AD libraries like ForwardDiff.jl and Zygote.jl is effortless and fast.","category":"page"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"using Transits\nusing Zygote\n\nlightcurve(X) = compute(PolynomialLimbDark(X[3:end]), X[1], X[2])\ngrad(X) = lightcurve'(X) # Zygote gradient\ngrad([0.1, 0.1, 0.4, 0.26])\n\n# output\n4-element Vector{Float64}:\n  0.0004972185834858653\n -0.2419262730830416\n -0.0048107583897073185\n -0.0024501564976671724","category":"page"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"To help demonstrate the logic behind these chain rules, here we derive a simple gradient function manually.","category":"page"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"using ChainRules\n\nu_n = [0.4, 0.26]\nμ = 0.1\nror = 0.1\nX0 = [μ, ror, u_n...]\n\nfunction gradr(X)\n    ld, ld_pullback = rrule(PolynomialLimbDark, X[3:end])\n    f, f_pullback = rrule(compute, ld, X[1], X[2])\n\n    f̄ = one(eltype(X))\n    _, l̄d, b̄, r̄ = f_pullback(f̄)\n    _, ū_n = ld_pullback(l̄d)\n    return [b̄, r̄, ū_n...]\nend\n\ngradr([0.1, 0.1, 0.4, 0.26])\n\n# output\n4-element Vector{Float64}:\n  0.0004972185834858653\n -0.2419262730830416\n -0.0048107583897073185\n -0.0024501564976671724","category":"page"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"For the most granular support for gradients and jacobians, peer into the depths of polynomial/poly-grad.jl and polynomial/quad-grad.jl. These functions are not part of the public API and are not guaranteed any stability according to semantic versioning.","category":"page"},{"location":"api/#Orbits","page":"API/Reference","title":"Orbits","text":"","category":"section"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"SimpleOrbit\nKeplerianOrbit\nOrbits.relative_position","category":"page"},{"location":"api/#Transits.Orbits.SimpleOrbit","page":"API/Reference","title":"Transits.Orbits.SimpleOrbit","text":"SimpleOrbit(; period, duration, t0=0, b=0.0)\n\nCircular orbit parameterized by the basic observables of a transiting system.\n\nParameters\n\nperiod - The orbital period of the planets, nominally in days\nduration - The duration of the transit, similar units as period.\nt0 - The midpoint time of the reference transit, similar units as period\nb - The impact parameter of the orbit, unitless\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.Orbits.KeplerianOrbit","page":"API/Reference","title":"Transits.Orbits.KeplerianOrbit","text":"KeplerianOrbit(; kwargs...)\n\nKeplerian orbit parameterized by the basic observables of a transiting 2-body system.\n\nParameters\n\nperiod/P – The orbital period of the planet [d].\nt0/t_0 – The midpoint time of the reference transit [d].\ntp/t_p – The time of periastron [d].\nduration/τ/T – The transit duration [d].\na – The semi-major axis [R⊙].\naR_star/aRs – The ratio of the semi-major axis to star radius.\nR_planet/Rp – The radius of the planet [R⊙].\nR_star/Rs – The radius of the star [R⊙].\nrho_star/ρ_star – The spherical star density [M⊙/R⊙³].\nr/RpRs – The ratio of the planet radius to star radius.\nb – The impact parameter, bounded between 0 ≤ b ≤ 1.\necc/e – The eccentricity of the closed orbit, bounded between 0 ≤ ecc < 1.\nincl – The inclination of the orbital plane relative to the axis perpendicular to the          reference plane [rad]\nomega/ω – The argument of periapsis [rad].\ncos_omega/cos_ω – The cosine of the argument of periapsis.\nsin_omega/sin_ω – The sine of the argument of periapsis.\nOmega/Ω – The longitude of the ascending node [rad].\nM_planet/Mp – The mass of the planet [M⊙].\nM_star/Ms – The mass of the star [M⊙].\n\nValid combinations\n\nThe following flowchart can be used to determine which parameters can define a KeplerianOrbit:\n\nThe period or a must be given. If both given, then neither M_star or rho_star can be defined because the stellar density is now implied.\nOnly incl or b can be given.\nIf ecc is given, then omega must also be given.\nIf no stellar parameters are given, the central body is assumed to be the Sun. If only rho_star is given, then R_star is defined to be 1 solar radius. Otherwise, at most two of M_star, R_star, and rho_star can be given.\nEither t0 or tp must be given, but not both.\n\n\n\n\n\n","category":"type"},{"location":"api/#Transits.Orbits.relative_position","page":"API/Reference","title":"Transits.Orbits.relative_position","text":"relative_position(::AbstractOrbit, t)\n\nThe relative position, [x, y, z], of the companion compared to the host at time t. In other words, this is the vector pointing from the host to the companion along the line of sight. Nominally, the units of this distance are relative to the host's radius. For example, a distance of 2 is 2 stellar radii.\n\n\n\n\n\n","category":"function"},{"location":"api/#Distributions","page":"API/Reference","title":"Distributions","text":"","category":"section"},{"location":"api/","page":"API/Reference","title":"API/Reference","text":"Kipping13","category":"page"},{"location":"api/#Transits.Kipping13","page":"API/Reference","title":"Transits.Kipping13","text":"Kipping13()\n\nA non-informative prior for two-parameter limb-darkening coefficients using triangular sampling (Kipping 2013).\n\nExamples\n\njulia> using Random; rng = Random.seed!(10);\n\njulia> rand(rng, Kipping13())\n2-element Vector{Float64}:\n 0.24716310305467298\n 0.08836997249298882\n\njulia> rand(rng, Kipping13(), 5)\n2×5 Matrix{Float64}:\n 0.0664907  0.124817   1.00732    0.806902  0.74165\n 0.520411   0.222718  -0.389412  -0.314755  0.0768429\n\nReferences\n\nKipping (2013)\"Efficient, uninformative sampling of limb darkening coefficients for two-parameter laws\"\n\n\n\n\n\n","category":"type"},{"location":"introduction/#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"introduction/#Historical-overview","page":"Introduction","title":"Historical overview","text":"","category":"section"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"Transit light curves are an essential tool used for the detection of exoplanets. To date, there have been over 4,300 confirmed planets discovered in over 3,400 different star systems, with an additional 2,400 candidates currently awaiting follow-up analysis and validation[1]. Since the first confirmed discovery of an exoplanet – as part of a multi-planetary system in 1992[2], and the first exoplanet discovered around a Sun-like star shortly after in 1995[3] – there has been an explosion in new discoveries, thanks in large part to the successful Kepler/K2 and TESS space missions. The large majority of these planets have been detected via the transit method:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: ) Exoplanet Archive","category":"page"},{"location":"introduction/#Transit-method","page":"Introduction","title":"Transit method","text":"","category":"section"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"This method works by observing the dimming in apparent brightness of a star as a planet passes in front of it from our point of view. The plot of the star's brightness as a function of time defines the white light curve as seen in the schematic below:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: )","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"\"How Do You Find an Exoplanet?\" by John Asher Johnson","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"Even just starting with a simple single planet system in a circular orbit, there is already a wealth of information encoded in this diagram. These observations give us insight not only into the bulk properties of the planet, but into the architecture of its orbital system and characteristics of its host star as well. For example, direct observables from the light curve like the transit duration (T) and ingress/egress time (tau) give us information about how tilted its orbit is and how fast the planet is traveling, while the transit depth (delta) gives us a direct measure of the size of the planet relative to its star. For circular orbits, these are nicely summarized by:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"beginaligned\nfracR_textpR_* = delta^12 \n\nb^2 = 1 - delta^12fracTt \n\nfracaR_* = fracPdelta^142pi\nleft(frac4Ttauright)^12 \n\nrho_* = frac3PGpi^2left(fracdelta^14sqrtTtauright)^3 quad\nendaligned","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"where P is the period of the planet's orbit and a its semi-major axis, b is the impact parameter, R_* is the radius of its star, and rho_* is the stellar density.","category":"page"},{"location":"introduction/#Limb-darkening","page":"Introduction","title":"Limb darkening","text":"","category":"section"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"Not shown above is an added dimension that Transits.jl excels in, limb darkening, demonstrated in the schematic below:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: )","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"ASTR 236 class notes","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"This effect is intimately related to the shape of the light curve, and allows us to constrain the brightness profile of the star itself. As we will see next, the method of transit light curves is not just useful for the detection of exoplanets, but also for taking it to the next step of characterizing its atmosphere.","category":"page"},{"location":"introduction/#Transmission-spectroscopy","page":"Introduction","title":"Transmission spectroscopy","text":"","category":"section"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"If we perform the technique of transit light curve modeling on a wavelength-by-wavelength basis, we can further probe the properties of the host star and begin to make predictions about the properties of the planet's atmosphere, such as its chemical composition and whether clouds/hazes are likely to be present at higher altitudes. This analysis begins in the same way as with the white light curve seen above, only now a wavelength binned light curve is measured at a range of different wavelengths:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: ) Adapted from Weaver et al. (2021, submitted)","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"Plotting these wavelength dependent transit depths then builds a transmission spectrum, which is filled with information about the planet's atmosphere and its star, summarized below:","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: Text here!) Benneke & Seager (2012)","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"(Image: ) Rackham, Apai, & Giampapa (2018)","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"Performing forward modeling (see, e.g., Kempton et al. 2016, Goyal et al. 2017) and retrievals (see, e.g., Barstow et al. 2020 and references therein) using these frameworks then allows us to explore exoplanetary atmospheres in never before seen detail.","category":"page"},{"location":"introduction/#Summary","page":"Introduction","title":"Summary","text":"","category":"section"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"The detection and characterization of exoplanets through their transit light curves is a relatively new technique in the field of astronomy, with recent advances only being made possible through novel uses of large, ground-based telescopes and soon in the future with planned ELTs and space based missions like JWST. Studies using these observing facilities will require the fast and precise computation of transit light curves, which Transits.jl aims to provide.","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"[1]: https://exoplanetarchive.ipac.caltech.edu/","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"[2]: https://ui.adsabs.harvard.edu/abs/1992Natur.355..145W/abstract","category":"page"},{"location":"introduction/","page":"Introduction","title":"Introduction","text":"[3]: https://ui.adsabs.harvard.edu/abs/1995Natur.378..355M/abstract","category":"page"},{"location":"bench/#Benchmarks","page":"Benchmarks","title":"Benchmarks","text":"","category":"section"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"Transits.jl aims to be at least as fast as similar tools. Limbdark.jl is also written in Julia and Agol et al. (2020) showed it outperforms starry, PyTransit, and batman in both runtime speed and numerical accuracy. The following benchmarks are works in progress, but they already show a marginal improvement on the Limbdark.jl implementation.","category":"page"},{"location":"bench/#Setup","page":"Benchmarks","title":"Setup","text":"","category":"section"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"warning: Warning\nThese benchmarks are works in progress","category":"page"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"The code can be found in bench/. You'll need to set up the environment yourself, including the installation of Limbdark.jl.","category":"page"},{"location":"bench/#Performance","page":"Benchmarks","title":"Performance","text":"","category":"section"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"(Image: )","category":"page"},{"location":"bench/#Comparison-with-Limbdark.jl","page":"Benchmarks","title":"Comparison with Limbdark.jl","text":"","category":"section"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"(Image: )","category":"page"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"","category":"page"},{"location":"bench/","page":"Benchmarks","title":"Benchmarks","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Transits","category":"page"},{"location":"#Transits.jl","page":"Home","title":"Transits.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: GitHub) (Image: Build Status) (Image: PkgEval) (Image: Coverage)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: License) (Image: DOI)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Transits.jl provides flexible and powerful occultation curves with limb darkening. The goals of this package are, in this order","category":"page"},{"location":"","page":"Home","title":"Home","text":"have a simple interface with high composability\nbe flexible with respect to numeric types and application\nbe fully compatible with ChainRules.jl automatic differentiation (AD) system to leverage the derived analytical gradients\nprovide a codebase that is well-organized, instructive, and easy to extend\nmaintain high performance: at least as fast as similar tools","category":"page"},{"location":"","page":"Home","title":"Home","text":"In particular, PolynomialLimbDark implements the \"starry\" limb darkening method, which solves the flux integral analytically. This provides floating-point errors and runtimes that are best in class.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install use Pkg. From the REPL, press ] to enter Pkg-mode","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Transits","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to use the most up-to-date version of the code, check it out from main","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Transits#main","category":"page"},{"location":"#Citations","page":"Home","title":"Citations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use Transits.jl or a derivative of it in your work please consider citing it at the Zenodo DOI. If you use PolynomialLimbDark or QuadLimbDark please also cite Agol et al. (2020) and Luger et al. (2019). If you use Kipping13 please cite Kipping (2013). BibTeX for all those citations can be found in CITATIONS.bib.","category":"page"}]
}
