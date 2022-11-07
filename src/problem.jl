struct Problem
    F::AbstractDatafit
    G::AbstractPenalty
    A::Matrix
    y::Vector
    λ::Float64
    a::Vector
    m::Int
    n::Int
    function Problem(
        F::AbstractDatafit,
        G::AbstractPenalty,
        A::Matrix,
        y::Vector,
        λ::Float64;
    )

        m, n = size(A)

        # Consistency checks
        (length(y) == m) || throw(ArgumentError("Shape of A and y missmatch"))
        (λ >= 0.0) || throw(ArgumentError("The parameter λ must be positive"))

        # Precompute the column norms of A
        a = [norm(ai, 2)^2 for ai in eachcol(A)]

        return new(F, G, A, y, λ, a, m, n)
    end
end

function Base.print(io::IO, problem::Problem)
    println(io, "L0-penalized problem")
    println(io, "  Datafit : $(problem.F)")
    println(io, "  Penalty : $(problem.G)")
    println(io, "  Dims    : $(problem.m) x $(problem.n)")
    println(io, "  λ       : $(round(problem.λ, digits=6)) ($(round(problem.λ / compute_λmax(problem.F, problem.G, problem.A, problem.y), digits=6))λmax)")
end

function objective(problem::Problem, x::Vector, Ax::Vector)
    FAx = value(problem.F, problem.y, Ax)
    Gx = value(problem.G, x)
    return FAx + problem.λ * (norm(x, 0) + Gx)
end

objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

function compute_λmax(F::AbstractDatafit, G::AbstractPenalty, A::Matrix, y::Vector)
    return norm(A' * gradient(F, y, zeros(length(y))), Inf) / G.τ
end