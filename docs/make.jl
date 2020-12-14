push!(LOAD_PATH,"../src/")

using Documenter
using HOODESolver

makedocs(
    sitename = "HOODESolver.jl",
    format = Documenter.HTML(),
    modules = [HOODESolver],
    pages = ["Documentation" => "index.md",
             "Types"         => "types.md",
             "Functions"     => "functions.md"],
    repo = "https://github.com/ymocquar/HOODESolver.jl/blob/{commit}{path}#{line}"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
