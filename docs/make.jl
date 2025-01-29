using NSInterpolate
using Documenter

DocMeta.setdocmeta!(NSInterpolate, :DocTestSetup, :(using NSInterpolate); recursive=true)

makedocs(;
    modules=[NSInterpolate],
    authors="Mark Dransfield",
    sitename="NSInterpolate.jl",
    format=Documenter.HTML(;
        canonical="https://bamburgh.github.io/NSInterpolate.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bamburgh/NSInterpolate.jl",
    devbranch="main",
)
