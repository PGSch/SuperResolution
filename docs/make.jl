using Documenter, SuperResolution

makedocs(;
    modules=[SuperResolution],
    sitename="SuperResolution.jl",
    authors="Patrick Schneider",
    #repo="https://github.com/PGSch/SuperResolution.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(prettyurls=false),#;
    #     prettyurls=get(ENV, "CI", "false") == "true",
    #     canonical="https://PGSch.github.io/SuperResolution.jl",
    #     assets=String[],
    # ),
    pages=[
        "Home" => "index.md",
       #  "Algorithms" => [
       #     "algorithms.md",
       #     "init.md",
       #     "kmeans.md",
       #     "kmedoids.md",
       #     "hclust.md",
       #     "mcl.md",
       #     "affprop.md",
       #     "dbscan.md",
       #     "fuzzycmeans.md",
       # ],
       # "validate.md",
    ],
)

deploydocs(;
    repo="github.com/PGSch/SuperResolution.jl",
)
