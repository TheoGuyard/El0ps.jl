push!(LOAD_PATH, "../src/")
using Documenter, El0ps

makedocs(
    sitename = "El0ps.jl",
    pages    = [
        "Home" => "index.md",
        "Manual" => [
            "Quick start" => "manual/quickstart.md",
            "Solving problems" => "manual/optimize.md",
            "Fitting paths" => "manual/path.md",
            "Custom functions" => "manual/custom.md",
        ],
        "Library" => [
            "Problem" => "library/problem.md",
            "Datafits" => "library/datafits.md",
            "Penalties" => "library/penalties.md",
            "Solvers" => "library/solvers.md",
            "Bounding step" => "library/bounding.md",
            "Regularization path" => "library/path.md",
        ]
    ]
)

deploydocs(
    repo = "github.com/TheoGuyard/El0ps.jl.git",
)
