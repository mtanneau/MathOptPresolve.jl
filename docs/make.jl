using MathOptPresolve
using Documenter

makedocs(;
    modules=[MathOptPresolve],
    authors="Mathieu Tanneau and contributors",
    repo="https://github.com/mtanneau/MathOptPresolve.jl/blob/{commit}{path}#L{line}",
    sitename="MathOptPresolve.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mtanneau.github.io/MathOptPresolve.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "User guide" => [
            "Rules" => "manual/rules.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/mtanneau/MathOptPresolve.jl",
)
