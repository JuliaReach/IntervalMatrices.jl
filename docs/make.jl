using Documenter, IntervalMatrices

DocMeta.setdocmeta!(IntervalMatrices, :DocTestSetup, :(using IntervalMatrices);
                    recursive=true)

makedocs(
    sitename = "IntervalMatrices.jl",
    modules = [IntervalMatrices],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/aligned.css"], sidebar_sitename=true),
    strict = true,
    pages = [
        "Home" => "index.md",
        "Library" => Any[
        "Types" => "lib/types.md",
        "Methods" => "lib/methods.md"],
        "About" => "about.md",
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/IntervalMatrices.jl.git",
    push_preview = true,
)
