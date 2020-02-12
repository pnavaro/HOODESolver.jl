push!(LOAD_PATH,"../src/")

using Documenter
using HighlyOscillatoryProblems

makedocs(
    sitename = "HighlyOscillatoryProblems.jl",
    format = Documenter.HTML(),
    modules = [HighlyOscillatoryProblems],
    pages = ["Documentation" => "index.md",
             "Types"         => "types.md",
             "Functions"     => "functions.md"],
    repo = "https://gitlab.inria.fr/ua/HighlyOscillatoryProblems.jl/blob/{commit}{path}#{line}"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
