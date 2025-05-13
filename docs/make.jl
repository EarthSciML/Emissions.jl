using Emissions
using Documenter

DocMeta.setdocmeta!(Emissions, :DocTestSetup, :(using Emissions); recursive = true)

makedocs(;
    modules = [Emissions],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/ctessum/Emissions.jl/blob/{commit}{path}#{line}",
    sitename = "Emissions.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ctessum.github.io/Emissions.jl",
        assets = String[]
    ),
    pages = ["Home" => "index.md"]
)

deploydocs(; repo = "github.com/ctessum/Emissions.jl", devbranch = "main")
