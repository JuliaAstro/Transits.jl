
# Getting Started


## Usage

```julia
using Transits

orbit = SimpleOrbit(period=3, duration=1)
u = [0.4, 0.26] # quad limb dark
ld = PolynomialLimbDark(u)

t = range(-1, 1, length=1000) # days from t0
rs = range(0, 0.2, length=10) # radius ratio

fluxes = @. ld(orbit, t, rs')
```

![](https://github.com/JuliaAstro/Transits.jl/raw/main/limbdark.png)

## Integrated and Secondary Curves

[`IntegratedLimbDark`](@ref) can be used to numerically integrate each light curve exposure in time

```julia
ld = IntegratedLimbDark([0.4, 0.26])
orbit = SimpleOrbit(period=3, duration=1)
t = range(-1, 1, length=1000)
texp = [0.1 0.2 0.3]
# no extra calculations made
flux = @. ld(orbit, t, 0.2)
# use quadrature to find time-averaged flux for each t
flux_int = @. ld(orbit, t, 0.2, texp) 
```

![](https://github.com/JuliaAstro/Transits.jl/raw/main/integrated.png)

[`SecondaryLimbDark`](@ref) can be used to generate secondary eclipses given a surface brightness ratio

```julia
ld = SecondaryLimbDark([0.4, 0.26], brightness_ratio=0.1)
ld_int = IntegratedLimbDark(ld) # composition works flawlessly

orbit = SimpleOrbit(period=4, duration=1)
t = range(-1.25, 2.75, length=1000)
rs = range(0.01, 0.1, length=6)

f = @. ld(orbit, t, rs')
f_int = @. ld_int(orbit, t, rs', texp=0.3)
```

![](https://github.com/JuliaAstro/Transits.jl/raw/main/secondary.png)

## Using Units

Units from `Unitful.jl` are a drop-in substitution for numbers

```julia
using Unitful
orbit = SimpleOrbit(period=10u"d", duration=5u"hr")
t = range(-6, 6, length=1000)u"hr"
flux = @. ld(orbit, t, 0.1)
```
