using IntervalObservers
using Documenter

DocMeta.setdocmeta!(IntervalObservers, :DocTestSetup, :(using IntervalObservers); recursive=true)

# Build to Hugo static folder for integration with the website
hugo_static_path = joinpath(@__DIR__, "..", "..", "IntervalObserversWebsite", "static", "api")

makedocs(;
    modules=[IntervalObservers],
    authors="ousmane-prog <ousmane-junior.sane@etu.univ-amu.fr>",
    sitename="IntervalObservers.jl API Documentation",
    format=Documenter.HTML(;
        canonical="https://yourdomain.com/api/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "API Reference" => "api.md",
    ],
    build=hugo_static_path,
)

deploydocs(;
    repo="github.com/ousmane-prog/IntervalObservers.jl",
    devbranch="main",
)
