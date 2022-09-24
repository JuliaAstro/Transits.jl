using Plots

plotsdir(args...) = joinpath(@__DIR__, "output", args...)
isdir(plotsdir()) && rm(plotsdir(); recursive=true)
mkpath(plotsdir())

function mnerror(prec_frac, prec_abs; kwargs...)
    term = 1:length(prec_frac)
    m1 = @. !iszero(prec_abs)
    prec_abs = prec_abs[m1]
    prec_frac = filter(!iszero, prec_frac)

    plot(
        term[m1],
        prec_abs;
        xlabel="term",
        ylabel="abs. err",
        label="",
        yscale=:log10,
        rightmargin=70Plots.px,
        kwargs...,
    )
    return plot!(
        twinx(),
        prec_frac;
        yscale=:log10,
        ylabel="frac. err",
        label="",
        ygrid=false,
        kwargs...,
    )
end
