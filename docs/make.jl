using Documenter
using DocumenterCitations
using Plots
using HOODESolver

ENV["GKSwstype"] = "100"

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))
makedocs(
    sitename = "HOODESolver.jl",
    authors="Yves Mocquard, Nicolas Crouseilles and Pierre Navaro",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pnavaro.github.io/HOODESolver.jl",
        assets=String[],
    ),
    modules = [HOODESolver],
    pages = ["Documentation"    => "index.md",
             "Numerical Method" => "numerical_method.md",
             "DifferentialEquations"  => "common_interface.md",
             "Quickstart"       => "quickstart.md",
             "Charged Particle" => "charged_particle.md",
             "Future work"      => "future_work.md",
             "Types"            => "types.md",
             "Functions"        => "functions.md"],
    repo = "https://github.com/pnavaro/HOODESolver.jl/blob/{commit}{path}#{line}",
    plugins = [bib], 
)

deploydocs(;
branch = "gh-pages",
    devbranch = "master",
    repo="github.com/pnavaro/HOODESolver.jl",
)
