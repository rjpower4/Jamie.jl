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
            "Utilities" => "util.md",
        ],
        "CRTBP" => "crtbp.md",
    ],
)

deploydocs(
    repo = "github.com/rjpower4/Jamie.jl.git",
)
