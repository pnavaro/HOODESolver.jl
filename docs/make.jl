using Documenter
using HighlyOscillatoryProblems

makedocs(
    sitename = "HighlyOscillatoryProblems",
    format = Documenter.HTML(),
    modules = [HighlyOscillatoryProblems]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
