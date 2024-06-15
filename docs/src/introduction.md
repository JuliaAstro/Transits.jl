# Introduction

## Historical overview

Transit light curves are an essential tool used for the detection of
[exoplanets](https://en.wikipedia.org/wiki/Exoplanet). To date, there have been over 4,300
confirmed planets discovered in over 3,400 different star systems, with an additional
2,400 candidates currently awaiting follow-up analysis and validation.[^1] Since the first
confirmed discovery of an exoplanet -- as part of a multi-planetary system in 1992 [WF1992](@cite),
and the first exoplanet discovered around a Sun-like star shortly after in 1995 [MQ1995](@cite) --
there has been an explosion in new discoveries, thanks in large part to the successful
[Kepler/K2](https://www.nasa.gov/mission_pages/kepler/main/index.html) and
[TESS](https://tess.mit.edu/) space missions. The large majority of these planets have
been detected via the [transit method](https://exoplanets.nasa.gov/faq/31/whats-a-transit/):

![](https://exoplanetarchive.ipac.caltech.edu/exoplanetplots/exo_dischist_cumulative_cb.png)
[Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu)

## Transit method

This method works by observing the dimming in apparent brightness of a star as a planet
passes in front of it from our point of view. The plot of the star's brightness as a
function of time defines the *white light curve* as seen in the schematic below:

![](https://upload.wikimedia.org/wikipedia/commons/1/10/Theoretical_Transiting_Exoplanet_Light_Curve.jpg)

["How Do You Find an Exoplanet?" by John Asher Johnson](https://www.google.com/books/edition/How_Do_You_Find_an_Exoplanet/-DNJCgAAQBAJ?hl=en)

Even just starting with a simple single planet system in a circular orbit, there is
already a wealth of information encoded in this diagram. These observations give us
insight not only into the bulk properties of the planet, but into the architecture of its
orbital system and characteristics of its host star as well. For example, *direct
observables* from the light curve like the *transit duration* $(T)$ and *ingress/egress*
time $(\tau)$ give us information about how tilted its orbit is and how fast the planet is
traveling, while the *transit depth* $(\delta)$ gives us a direct measure of the size of
the planet relative to its star. For circular orbits, these are nicely summarized by:

```math
\begin{aligned}
\frac{R_\text{p}}{R_*} &= \delta^{1/2} \\

b^2 &= 1 - \delta^{1/2}\frac{T}{t} \\

\frac{a}{R_*} &= \frac{P\delta^{1/4}}{2\pi}
\left(\frac{4}{T\tau}\right)^{1/2} \\

\rho_* &= \frac{3P}{G\pi^2}\left(\frac{\delta^{1/4}}{\sqrt{T\tau}}\right)^3 \quad,
\end{aligned}
```

where $P$ is the period of the planet's orbit and $a$ its semi-major axis, $b$ is the
impact parameter, $R_*$ is the radius of its star, and $\rho_*$ is the stellar
density.

### Limb darkening
Not shown above is an added dimension that `Transits.jl` excels in,
[limb darkening](https://en.wikipedia.org/wiki/Limb_darkening),
demonstrated in the schematic below:

![](https://user-images.githubusercontent.com/25312320/108404912-712f1c00-71ee-11eb-968e-b34001fe7a55.jpg)

[ASTR 236 class notes](http://www.astro.utoronto.ca/~astrolab/files/AST326_LimbDarkening_2017.pdf)

This effect is intimately related to the shape of the light curve, and allows us to
constrain the brightness profile of the star itself. As we will see next, the method of
transit light curves is not just useful for the detection of exoplanets, but also for
taking it to the next step of characterizing its atmosphere.

## Transmission spectroscopy
If we perform the technique of transit light curve modeling on a wavelength-by-wavelength
basis, we can further probe the properties of the host star and begin to make predictions
about the properties of the planet's atmosphere, such as its chemical composition and
whether clouds/hazes are likely to be present at higher altitudes. This analysis begins in
the same way as with the white light curve seen above, only now a *wavelength binned light
curve* is measured at a range of different wavelengths:

![](https://user-images.githubusercontent.com/25312320/108020235-f1386480-6fe9-11eb-87f2-4970dabd7839.png)
Adapted from Weaver et al. (2021, *submitted*)

Plotting these wavelength dependent transit depths then builds a *transmission
spectrum*, which is filled with information about the planet's atmosphere and its star,
summarized below:

![Text here!](https://user-images.githubusercontent.com/25312320/108021680-124e8480-6fed-11eb-8eaf-bbf9b0df217b.jpg)
[BS2012](@citet)

![](https://s3.amazonaws.com/aasie/images/0004-637X/853/2/122/apjaaa08cf1_hr.jpg)
[RAG2018](@citet)

Performing forward modeling (see, e.g., [KLOASC2017](@citet),
[GMSDTAECSBNMCH2018](@citet)) and retrievals (see, e.g., [BCGLRW2020](@citet)
and references therein) using these frameworks then allows us to explore
exoplanetary atmospheres in never before seen detail.

## Summary
The detection and characterization of exoplanets through their transit light curves is a relatively new technique
in the field of astronomy, with recent advances only being made possible through novel
uses of large, ground-based telescopes and soon in the future with planned
[ELTs](https://en.wikipedia.org/wiki/Extremely_large_telescope) and space based
missions like [JWST](https://www.jwst.nasa.gov/). Studies using these observing facilities
will require the fast and precise computation of transit light curves, which
[`Transits.jl`](https://github.com/JuliaAstro/Transits.jl)
aims to provide.

## Footnotes and references
[^1]: https://exoplanetarchive.ipac.caltech.edu/

```@bibliography
Pages = ["introduction.md"]
```
