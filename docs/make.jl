push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using Plots
using HOODESolver

ENV["GKSwstype"] = "100"

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))
makedocs(
    bib, 
    sitename = "HOODESolver.jl",
    authors="Yves Mocquard",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ymocquar.github.io/HOODESolver.jl",
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
    repo = "https://github.com/ymocquar/HOODESolver.jl/blob/{commit}{path}#{line}"
)

deploydocs(;
branch = "gh-pages",
    devbranch = "master",
    repo="github.com/ymocquar/HOODESolver.jl",
)
