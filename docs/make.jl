using Documenter
using HighlyOscillatoryProblems

makedocs(
    sitename = "HighlyOscillatoryProblems.jl",
    format = Documenter.HTML(),
    modules = [HighlyOscillatoryProblems],
    pages = ["Documentation" => "index.md",
             "Types"         => "types.md",
             "Functions"     => "functions.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
