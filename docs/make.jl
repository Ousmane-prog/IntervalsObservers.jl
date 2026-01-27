using IntervalsObservers
using Documenter

DocMeta.setdocmeta!(IntervalsObservers, :DocTestSetup, :(using IntervalsObservers); recursive=true)

makedocs(;
    modules=[IntervalsObservers],
    authors="ousmane-prog <ousmane-junior.sane@etu.univ-amu.fr>",
    sitename="IntervalsObservers.jl",
    format=Documenter.HTML(;
        canonical="https://ousmane-prog.github.io/IntervalsObservers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ousmane-prog/IntervalsObservers.jl",
    devbranch="main",
)
