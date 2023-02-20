push!(LOAD_PATH, "../src/")
using Documenter, El0ps

makedocs(
    sitename = "El0ps.jl",
    pages    = [
        "Home" => "index.md",
        "Manual" => [
            "Quick start" => "manual/quickstart.md",
            "Creating problems" => "manual/problem.md",
            "Solving problems" => "manual/optimize.md",
            "Fitting paths" => "manual/path.md",
            "Custom functions" => "manual/custom.md",
        ],
        "Library" => [
            "Problem" => "library/problem.md",
            "Datafits" => "library/datafits.md",
            "Perturbations" => "library/perturbations.md",
            "Solver" => "library/solver.md",
            "Path" => "library/path.md",
        ]
    ]
)

deploydocs(
    repo = "github.com/TheoGuyard/El0ps.jl.git",
)
