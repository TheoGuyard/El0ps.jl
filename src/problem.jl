"""
    Problem

Concrete type for a problem of the form 
    
    min F(y,Ax) + λ * (||x||_0 + G(x))

where A ∈ R^{m×n}, y ∈ R^m and λ > 0.
"""
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

"""
    objective(problem::Problem, x::Vector, Ax::Vector)

Computes the value of F(y,Ax) + λ * (||x||_0 + G(x)) with Ax already given.
"""
function objective(problem::Problem, x::Vector, Ax::Vector)
    FAx = value(problem.F, problem.y, Ax)
    Gx = value(problem.G, x)
    return FAx + problem.λ * (norm(x, 0) + Gx)
end

"""
    objective(problem::Problem, x::Vector)

Computes the value of F(y,Ax) + λ * (||x||_0 + G(x)).
"""
objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

"""
    compute_λmax(F::AbstractDatafit, G::AbstractPenalty, A::Matrix, y::Vector)

Return a value of λ such that 0 is a solution of min F(y,Ax) + λ * (||x||_0 + G(x)).
"""
function compute_λmax(F::AbstractDatafit, G::AbstractPenalty, A::Matrix, y::Vector)
    return norm(A' * gradient(F, y, zeros(length(y))), Inf) / G.τ
end