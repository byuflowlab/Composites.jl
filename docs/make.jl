using Documenter, Composites

makedocs(
    sitename = "Composites",
    modules = [Composites],
    pages   = ["Home" => "index.md",
        "Library" => "library.md"]
)

deploydocs(
    repo   = "github.com/byuflowlab/Composites.jl")
