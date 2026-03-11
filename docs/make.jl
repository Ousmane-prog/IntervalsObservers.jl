using IntervalObservers
using Documenter

# For reproducibility
mkpath(joinpath(@__DIR__, "src", "assets"))
cp(
    joinpath(@__DIR__, "Manifest.toml"),
    joinpath(@__DIR__, "src", "assets", "Manifest.toml");
    force=true,
)
cp(
    joinpath(@__DIR__, "Project.toml"),
    joinpath(@__DIR__, "src", "assets", "Project.toml");
    force=true,
)

repo = "github.com/ousmane-prog/IntervalObservers.jl.git"
repo_link = "https://github.com/ousmane-prog/IntervalObservers.jl"

makedocs(
    modules = [IntervalObservers],
    authors = "ousmane-prog <ousmane-junior.sane@etu.univ-amu.fr>",
    sitename = "IntervalObservers.jl",
    repo = repo_link,
    format = Documenter.HTML(
        repolink = repo_link,
        edit_link = "main",
        prettyurls = true
    ),
    doctest = false,
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Nonlinear System Example" => "nonlinear_system_example.md",
        # "API Reference" => "api.md",
    ],
)

deploydocs(
    repo = repo,
    devbranch = "main",
)