using Documenter, IntervalMatrices

makedocs(
    doctest = true,  # use this flag to skip doctests (saves time!)
    modules = [IntervalMatrices],
    format = :html,
    assets = ["assets/juliareach.css"],
    sitename = "IntervalMatrices.jl",
    pages = [
        "Home" => "index.md",
        "Library" => Any[
        "Types" => "lib/types.md",
        "Methods" => "lib/methods.md"],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/IntervalMatrices.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = nothing,
    make = nothing
)
