using Documenter
using Orbits
using Transits
using Documenter.Remotes: GitHub

setup = quote
    using Transits
end

DocMeta.setdocmeta!(Transits, :DocTestSetup, setup; recursive=true)

# gives us `pages` and `bib`
include("pages.jl")

makedocs(;
    modules=[Transits, Orbits],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo=GitHub("JuliaAstro/Transits.jl"),
    sitename="Transits.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.org/Transits/stable/",
        assets=String[],
    ),
    pages=pages,
    warnonly=[:missing_docs],
    plugins=[bib],
)

deploydocs(;
    repo = "github.com/JuliaAstro/Transits.jl",
    push_preview = true,
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
