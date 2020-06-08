push!(LOAD_PATH, "../src/")

using Documenter
using Jamie

makedocs(
    sitename="Jamie Documentation",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Overview" => [
            "Introduction" => "index.md",
            "Position-Velocity" => "posvel.md",
            "Celestial Bodies" => "body.md",
        ],
        "CRTBP" => "crtbp.md",
    ],
)
