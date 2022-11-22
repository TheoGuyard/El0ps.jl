push!(LOAD_PATH,"../src/")
using Documenter, El0ps

makedocs(
    sitename = "El0ps.jl",
    pages    = [
        "Home" => "index.md",
        "Examples" => [
            "Solve a problem" => "examples/optimize.md",
            "Fit a regularization path" => "examples/path.md",
        ],
        "Library" => [
            "Data generation" => "library/data.md",
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
