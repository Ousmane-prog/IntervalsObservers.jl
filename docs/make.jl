using IntervalObservers
using Documenter
using DocumenterCitations


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

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)
repo = "github.com/Ousmane-prog/IntervalsObservers.jl.git"
repo_link = "https://github.com/Ousmane-prog/IntervalsObservers.jl"

makedocs(
    authors = "Ousmane-prog <myboxe2@gmail.com>",
    sitename = "IntervalsObservers.jl",
    repo = repo_link,
    format = Documenter.HTML(
        repolink = repo_link,
        edit_link = "main",
        example_size_threshold = 65536,
        assets = String["../assets/custom.css"],
    ),
    plugins = [bib],
    doctest = false,
    draft = false,
    checkdocs = :exports,
    warnonly = [:autodocs_block],
    pages = [
        "Getting Started" => "getting_started.md",
        "Theory" => "theory.md",
        "Applications" => [
            "Case Study 1" => "applications.md",
            "Case Study 2" => "applications2.md",
        ],
        "Examples" => [
            "Nonlinear System Example" => "nonlinear_system_example.md",
        ],
        "API Reference" => "api.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = repo,
        devbranch = "main",
    )
end