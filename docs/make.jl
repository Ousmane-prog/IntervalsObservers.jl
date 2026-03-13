using IntervalObservers
using Documenter
using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const GENERATED_DIR = joinpath(@__DIR__, "src", "generated")

mkpath(GENERATED_DIR)

# Convert every Literate example to markdown
example_files = sort(filter(f -> endswith(f, ".jl"), readdir(EXAMPLES_DIR)))

generated_pages = Pair{String,String}[]

for file in example_files
    input_path = joinpath(EXAMPLES_DIR, file)

    # Generate markdown into docs/src/generated
    Literate.markdown(input_path, GENERATED_DIR; execute=true)

    md_file = replace(file, ".jl" => ".md")

    # Turn filename into a readable sidebar title
    title = replace(replace(file, ".jl" => ""), "_" => " ")
    title = titlecase(title)

    push!(generated_pages, title => joinpath("generated", md_file))
end

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
    pages = [
        "Home" => "index.md",
        "Getting Started" => "Getting_started.md",
        "Examples" => generated_pages,
    ],
)

deploydocs(
    repo = repo,
    devbranch = "main",
)