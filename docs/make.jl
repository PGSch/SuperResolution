using SuperResolution
using Documenter

makedocs(;
    modules=[SuperResolution],
    authors="Patrick Schneider",
    repo="https://github.com/PGSch/SuperResolution.jl/blob/{commit}{path}#L{line}",
    sitename="SuperResolution.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://PGSch.github.io/SuperResolution.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PGSch/SuperResolution.jl",
)
