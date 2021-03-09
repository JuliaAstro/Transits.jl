using Transits
using Documenter

setup = quote
    using Transits
end

DocMeta.setdocmeta!(Transits, :DocTestSetup, setup; recursive = true)

makedocs(;
    modules=[Transits],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/juliaastro/Transits.jl/blob/{commit}{path}#L{line}",
    sitename="Transits.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/Transits.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Getting Started" => "gettingstarted.md",
        "Benchmarks" => "bench.md",
        "API/Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/Transits.jl",
    push_preview=true,
    devbranch="main"
)
