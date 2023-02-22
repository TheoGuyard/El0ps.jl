"""
    Problem

Structure for L0-penalized problems.
"""
struct Problem
    f::AbstractDatafit
    h::AbstractPerturbation
    A::Matrix
    λ::Float64
    a::Vector
    m::Int
    n::Int
    λmax::Float64
end

"""
    Problem(
        f::AbstractDatafit,
        h::AbstractPerturbation,
        A::Matrix,
        λ::Float64,
    )

Instantiate a [`Problem`](@ref) of the form `min f(Ax) + λ (norm(x,0) + h(x))`.
"""
function Problem(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix, λ::Float64)
    m = size(A, 1)
    n = size(A, 2)
    a = [norm(ai, 2)^2 for ai in eachcol(A)]
    @assert dim_input(f) == m
    @assert λ >= 0.0
    @assert !any(a .≈ 0.0)
    λmax = compute_λmax(f, h, A)
    return Problem(f, h, A, λ, a, m, n, λmax)
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

Value of the objective of a [`Problem`](@ref) when `Ax` is already computed.
"""
function objective(problem::Problem, x::Vector, Ax::Vector)
    fval = value(problem.f, Ax)
    hval = value(problem.h, x)
    return fval + problem.λ * (norm(x, 0) + hval)
end

"""
    objective(problem::Problem, x::Vector)

Value of the objective of a [`Problem`](@ref).
"""
objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

"""
    compute_λmax(
        f::AbstractDatafit, 
        h::AbstractPerturbation, 
        A::Matrix, 
        y::Vector
    )

Return a value of `λ` such that `0` is a solution of a [`Problem`](@ref).
"""
function compute_λmax(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix)
    return norm(A' * gradient(f, zeros(dim_input(f))), Inf) / h.τ
end
