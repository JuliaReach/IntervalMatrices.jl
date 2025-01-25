using Documenter, IntervalMatrices, DocumenterCitations

DocMeta.setdocmeta!(IntervalMatrices, :DocTestSetup, :(using IntervalMatrices);
                    recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(; sitename="IntervalMatrices.jl",
         modules=[IntervalMatrices],
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true",
                                assets=["assets/aligned.css", "assets/citations.css"],
                                sidebar_sitename=true),
         pagesonly=true,
         plugins=[bib],
         pages=["Home" => "index.md",
                "Library" => Any["Types" => "lib/types.md",
                                 "Methods" => "lib/methods.md"],
                "About" => "about.md",
                "Bibliography" => "bibliography.md"])

deploydocs(; repo="github.com/JuliaReach/IntervalMatrices.jl.git",
           push_preview=true)
