using Emissions
using Documenter

DocMeta.setdocmeta!(Emissions, :DocTestSetup, :(using Emissions); recursive = true)

makedocs(;
    modules = [Emissions],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/EarthSciML/Emissions.jl/blob/{commit}{path}#{line}",
    sitename = "Emissions.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://emissions.earthsci.dev",
        assets = String[],
        repolink = "https://github.com/EarthSciML/Emissions.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Complete Tutorial" => "tutorial.md",
        "NEI Processing" => "nei_processing.md",
    ],
    checkdocs = :none
)

deploydocs(; repo = "github.com/EarthSciML/Emissions.jl.git")
