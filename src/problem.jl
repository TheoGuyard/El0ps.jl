"""
    Problem

Structure for L0-penalized problems.
"""
struct Problem
    f::AbstractDatafit
    h::AbstractPerturbation
    A::Matrix
    λ::Float64
    τ::Float64
    μ::Float64
    a::Vector
    m::Int
    n::Int
    λmax::Float64
end

"""
    Problem(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix, λ::Float64)

Instantiate a [`Problem`](@ref) of the form `min f(Ax) + λ * norm(x,0) + h(x)`.
"""
function Problem(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix, λ::Float64)
    m = size(A, 1)
    n = size(A, 2)
    a = [norm(ai, 2)^2 for ai in eachcol(A)]
    @assert dim_input(f) == m
    @assert λ >= 0.0
    @assert !any(a .≈ 0.0)
    τ = compute_τ(h, λ)
    μ = compute_μ(h, λ)
    λmax = compute_λmax(f, h, A)
    return Problem(f, h, A, λ, τ, μ, a, m, n, λmax)
end

function Base.show(io::IO, problem::Problem)
    println(io, "L0-penalized problem")
    println(io, "  Datafit : $(problem.f)")
    println(io, "  Perturb : $(problem.h)")
    println(io, "  Dims    : $(problem.m) x $(problem.n)")
    println(io, "  λ       : $(round(problem.λ, digits=4))")
    @printf "  λ/λmax  : %.1e" problem.λ / problem.λmax
end

"""
    objective(problem::Problem, x::Vector, Ax::Vector)

Value of the objective of a [`Problem`](@ref) when `Ax` is already computed.
"""
function objective(problem::Problem, x::Vector, Ax::Vector)
    fval = value(problem.f, Ax)
    hval = value(problem.h, x)
    return fval + problem.λ * norm(x, 0) + hval
end

"""
    objective(problem::Problem, x::Vector)

Value of the objective of a [`Problem`](@ref).
"""
objective(problem::Problem, x::Vector) = objective(problem, x, problem.A * x)

"""
    compute_λmax(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix)

Return a value of `λ` such that the all-zero vector is a solution of a [`Problem`](@ref).
"""
function compute_λmax(f::AbstractDatafit, h::AbstractPerturbation, A::Matrix; ϵ::Float64=1e-8)

    g = norm(A' * gradient(f, zeros(dim_input(f))), Inf)

    (compute_τ(h, 0.) >= g) && return 0.

    # Find λa and λb such that (τ(λa) - g) * (τ(λb) - g) < 0
    λa = 1.
    τa = compute_τ(h, λa)
    λf = (τa <= g) ? 10. : 0.1
    λb = λf * λa
    τb = compute_τ(h, λb)
    while (τa - g) * (τb - g) >= 0.0
        λb = λf * λb
        τb = compute_τ(h, λb)
        (λb == Inf) && return Inf
    end

    # Bisection method on to find λmax in [λa,λb] such that τ(λmax) ≈ g
    while abs(λa - λb) > ϵ
        λc = (λa + λb) / 2.0
        τc = compute_τ(h, λc)
        if (τa - g) * (τc - g) < 0.0
            λb = λc
            τb = compute_τ(h, λb)
        else
            λa = λc
            τa = compute_τ(h, λa)
        end
    end

    λmax = (λa + λb) / 2.0

    return λmax
end
