using IntervalObservers
using Documenter

DocMeta.setdocmeta!(IntervalObservers, :DocTestSetup, :(using IntervalObservers); recursive=true)

makedocs(;
    modules=[IntervalObservers],
    authors="ousmane-prog <ousmane-junior.sane@etu.univ-amu.fr>",
    sitename="IntervalObservers.jl",
    format=Documenter.HTML(;
        canonical="https://ousmane-prog.github.io/IntervalObservers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ousmane-prog/IntervalObservers.jl",
    devbranch="main",
)
