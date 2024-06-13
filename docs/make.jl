using Documenter
using Orbits
using Transits
using Documenter.Remotes: GitHub
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
    repo=GitHub("JuliaAstro/Transits.jl"),
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
