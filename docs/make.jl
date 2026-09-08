CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using QuantumInputOutput
using SecondQuantizedAlgebra
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterInterLinks
using DocumenterLandingPage

using Plots
default(; fmt = :png)

ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

include("make_md_examples.jl")
include("pages.jl")

bib = CitationBibliography("src/refs.bib"; style = :authoryear)
links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "Documenter" => "https://documenter.juliadocs.org/stable/",
    "SecondQuantizedAlgebra" => "https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/objects.inv",
)

DocMeta.setdocmeta!(
    QuantumInputOutput,
    :DocTestSetup,
    :(using QuantumInputOutput, SecondQuantizedAlgebra);
    recursive = true,
)

makedocs(
    sitename = "QuantumInputOutput.jl",
    authors = "Orjan Ameye and contributors",
    modules = [QuantumInputOutput],
    format = Documenter.HTML(;
        canonical = "https://qojulia.github.io/QuantumInputOutput.jl",
        assets = [asset("assets/favicon.png", class = :ico, islocal = true)],
    ),
    pages = pages,
    plugins = [bib, CodeBlocks(), LandingPage(), links],
    clean = true,
    linkcheck = false,
    warnonly = :missing_docs,
    draft = false,
    doctest = false, # run in the test suite
    checkdocs = :exports,
)

if CI
    deploydocs(
        repo = "github.com/qojulia/QuantumInputOutput.jl",
        devbranch = "main",
        target = "build",
        branch = "gh-pages",
        push_preview = true,
        versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
    )
end
