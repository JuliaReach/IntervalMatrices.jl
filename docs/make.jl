using Documenter, IntervalMatrices

DocMeta.setdocmeta!(IntervalMatrices, :DocTestSetup, :(using IntervalMatrices);
                    recursive=true)

makedocs(
    sitename = "IntervalMatrices.jl",
    modules = [IntervalMatrices],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/juliareach.css"]),
    pages = [
        "Home" => "index.md",
        "Library" => Any[
        "Types" => "lib/types.md",
        "Methods" => "lib/methods.md"],
        "About" => "about.md"
    ],
    strict = true
)

deploydocs(
    repo = "github.com/JuliaReach/IntervalMatrices.jl.git"
)
