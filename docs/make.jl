using IntervalObservers
using Documenter


mkpath(joinpath(@__DIR__, "src", "assets"))

cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force = true,
)

cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force = true,
)

repo = "github.com/ousmane-prog/IntervalObservers.jl.git"
repo_link = "https://github.com/ousmane-prog/IntervalObservers.jl"

makedocs(
    authors = "ousmane-prog <ousmane-junior.sane@etu.univ-amu.fr>",
    sitename = "IntervalObservers.jl",
    repo = repo_link,
    format = Documenter.HTML(
        repolink = repo_link,
        edit_link = "main",
        assets = String[],
    ),
    doctest = false,
    draft = false,
    pages = [
        "Getting Started" => "getting_started.md",
        "Theory" => "theory.md",
        "Applications" => "applications.md",
        "Examples" => [
            "Nonlinear System Example" => "nonlinear_system_example.md",
        ],
    ],
)

deploydocs(
    repo = repo,
    devbranch = "main",
)