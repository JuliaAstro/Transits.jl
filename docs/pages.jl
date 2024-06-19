using DocumenterCitations

pages = [
    "Home" => "index.md",
    "Introduction" => "introduction.md",
    "Getting Started" => "gettingstarted.md",
    "Benchmarks" => "bench.md",
    "API/Reference" => "api.md",
]

# put bib here for juliaastro.github.io so that the builder for that site can
# pick it up
bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)
