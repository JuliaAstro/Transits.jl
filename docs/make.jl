using Documenter
using Orbits
using Transits
using DocumenterCitations

setup = quote
    using Transits
end

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

DocMeta.setdocmeta!(Transits, :DocTestSetup, setup; recursive=true)

include("pages.jl")

makedocs(;
    modules=[Transits, Orbits],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/juliaastro/Transits.jl/blob/{commit}{path}#L{line}",
    sitename="Transits.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/Transits.jl",
        assets=String[],
    ),
    pages=pages,
    warnonly=[:missing_docs],
    plugins=[bib],
)

deploydocs(; repo="github.com/JuliaAstro/Transits.jl", push_preview=true, devbranch="main")
