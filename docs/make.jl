using Documenter
using HighOscillatoryProblems

makedocs(
    sitename = "HighOscillatoryProblems",
    format = Documenter.HTML(),
    modules = [HighOscillatoryProblems]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
