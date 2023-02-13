"""
    Problem

Problem of the form 
    
``\\min f(\\mathbf{y},\\mathbf{Ax}) + \\lambda (\\|x\\|_0 + h(\\mathbf{x}))``

where ``\\mathbf{A} \\in \\mathbf{R}^{m \\times n}``, 
``\\mathbf{y} \\in \\mathbf{R}^{m}`` and ``\\lambda > 0``.
"""
struct Problem
    f::AbstractDatafit
    h::AbstractPerturbation
    A::Matrix
    y::Vector
    λ::Float64
    a::Vector
    m::Int
    n::Int
    λmax::Float64
    function Problem(
        f::AbstractDatafit,
        h::AbstractPerturbation,
        A::Matrix,
        y::Vector,
        λ::Float64,
    )
        m = size(A, 1)
        n = size(A, 2)
        a = [norm(ai, 2)^2 for ai in eachcol(A)]
        @assert length(y) == m
        @assert λ >= 0.
        @assert !any(a .≈ 0.)
        λmax = compute_λmax(f, h, A, y)
        return new(f, h, A, y, λ, a, m, n, λmax)
    end
end

function Base.show(io::IO, problem::Problem)
    println(io, "L0-penalized problem")
    println(io, "  Datafit : $(problem.f)")
    println(io, "  Perturb : $(problem.h)")
    println(io, "  Dims    : $(problem.m) x $(problem.n)")
    println(io, "  λ       : $(round(problem.λ, digits=4))")
    print(io, "  λ/λmax  : $(round(problem.λ/problem.λmax, digits=4))")
end

"""
    objective(problem::Problem, x::Vector, Ax::Vector)

Returns the value of the objective of a [`Problem`](@ref) when `Ax` is already
computed.
"""
function objective(problem::Problem, x::Vector, Ax::Vector)
    FAx = value(problem.f, problem.y, Ax)
    hx = value(problem.h, x)
    return FAx + problem.λ * (norm(x, 0) + hx)
end

"""
    objective(problem::Problem, x::Vector)

Returns the value of the objective of a [`Problem`](@ref).
"""
objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

"""
    compute_λmax(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix, y::Vector)

Return a value of `λ` such that `0` is a solution of a [`Problem`](@ref).
"""
function compute_λmax(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix, y::Vector)
    return norm(A' * gradient(f, y, zeros(length(y))), Inf) / h.τ
end
