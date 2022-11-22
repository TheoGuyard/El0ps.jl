"""
    Problem

Problem of the form 
    
``\\min F(\\mathbf{y},\\mathbf{Ax}) + \\lambda (\\|x\\|_0 + G(\\mathbf{x}))``

where ``\\mathbf{A} \\in \\mathbf{R}^{m \\times n}``, 
``\\mathbf{y} \\in \\mathbf{R}^{m}`` and ``\\lambda > 0``.
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
        λ::Float64,
    )

        m, n = size(A)
        @assert length(y) == m
        @assert λ >= 0.0
        a = [norm(ai, 2)^2 for ai in eachcol(A)]
        @assert !any(a .≈ 0.)

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

Returns the value of the objective of a [`Problem`](@ref) when `Ax` is already
computed.
"""
function objective(problem::Problem, x::Vector, Ax::Vector)
    FAx = value(problem.F, problem.y, Ax)
    Gx = value(problem.G, x)
    return FAx + problem.λ * (norm(x, 0) + Gx)
end

"""
    objective(problem::Problem, x::Vector)

Returns the value of the objective of a [`Problem`](@ref).
"""
objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

"""
    compute_λmax(F::AbstractDatafit, G::AbstractPenalty, A::Matrix, y::Vector)

Return a value of `λ` such that `0` is a solution of a [`Problem`](@ref).
"""
function compute_λmax(F::AbstractDatafit, G::AbstractPenalty, A::Matrix, y::Vector)
    return norm(A' * gradient(F, y, zeros(length(y))), Inf) / G.τ
end