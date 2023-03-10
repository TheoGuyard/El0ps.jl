"""
    Path

Regularization path of a [`Problem`](@ref), i.e., solutions for different
values of `λ`.
"""
Base.@kwdef mutable struct Path
    λ::Vector = Vector()
    λratio::Vector = Vector()
    x::Vector{Vector} = Vector()
    converged::Vector{Bool} = Vector{Bool}()
    solve_time::Vector = Vector()
    node_count::Vector{Int} = Vector{Int}()
    objective_value::Vector = Vector()
    datafit_value::Vector = Vector()
    penalty_value::Vector = Vector()
    support_size::Vector{Int} = Vector{Int}()
    cv_mean::Vector{Float64} = Vector{Float64}()
    cv_std::Vector{Float64} = Vector{Float64}()
end

"""
    PathOptions

Options for a [`Path`](@ref). The path is computed over a
logarithmically-spaced grid `λ ∈ [λratio_max, λratio_min] * λmax` of
`λratio_num` different values. The value of `λmax` is computed using
[`compute_λmax`](@ref).

# Arguments

- `λratio_max::Float64` : Maximum value of `λ/λmax`.
- `λratio_min::Float64` : Minimum value of `λ/λmax`.
- `λratio_num::Int` : Number of values of `λ` in the regularization path.
- `max_support_size::Int` : Stop the path fitting when a solution with support
size `max_support_size` is recovered.
- `stop_if_unsolved::Bool` : If `true`, stop the path fitting if the
[`Problem`](@ref) at some `λ` is unsolved.
- `compute_cv::Bool` : If `true`, compute the cross-validation error over
`nb_folds` folds for each solution obtained in the [`Path`](@ref).
- `nb_folds::Int` : Number of cross-validation folds.
- `verbosity::Bool` : Whether to displays outputs.
"""
struct PathOptions
    λratio_max::Float64
    λratio_min::Float64
    λratio_num::Int
    max_support_size::Int
    stop_if_unsolved::Bool
    compute_cv::Bool
    nb_folds::Int
    verbosity::Bool
    function PathOptions(;
        λratio_max::Float64 = 1.0,
        λratio_min::Float64 = 1e-2,
        λratio_num::Int = 21,
        max_support_size::Int = typemax(Int),
        stop_if_unsolved::Bool = true,
        compute_cv::Bool = true,
        nb_folds::Int = 10,
        verbosity::Bool = true,
    )

        @assert 0.0 <= λratio_min <= λratio_max <= 1.0
        @assert 0 <= λratio_num
        @assert 0 <= max_support_size
        @assert 0 <= nb_folds

        return new(
            λratio_max,
            λratio_min,
            λratio_num,
            max_support_size,
            stop_if_unsolved,
            compute_cv,
            nb_folds,
            verbosity,
        )
    end
end

const PATH_HEAD_STRING = " λ/λmax   Conv     Time     Fval     Hval   Nnz  CV mean ±  CV std"

function display_path_head()
    println(repeat("-", length(PATH_HEAD_STRING)))
    println(PATH_HEAD_STRING)
    println(repeat("-", length(PATH_HEAD_STRING)))
    return nothing
end

function display_path_tail()
    println(repeat("-", length(PATH_HEAD_STRING)))
end

function display_path_info(path::Path, i::Union{Int,Nothing} = nothing)
    i = isa(i, Nothing) ? length(path.λratio) : i
    @printf "%.1e" path.λratio[i]
    @printf "   %s" string(path.converged[i])
    @printf "  %7.2f" path.solve_time[i]
    @printf "  %7.2f" path.datafit_value[i]
    @printf "  %7.2f" path.penalty_value[i]
    @printf "  %4d" path.support_size[i]
    if isnan(path.cv_mean[i])
        print("      NaN ±     NaN")
    else
        @printf "  %.1e ± %.1e" path.cv_mean[i] path.cv_std[i]
    end
    println()
end

function fill_path!(
    path::Path,
    f::AbstractDatafit,
    h::AbstractPenalty,
    A::Matrix,
    λ::Float64,
    λratio::Float64,
    result::AbstractResult,
    options::PathOptions,
)
    if options.compute_cv
        cv_mean, cv_std = compute_cv_statistics(result.x, f, A, options.nb_folds)
    else
        cv_mean, cv_std = NaN, NaN
    end
    push!(path.λ, λ)
    push!(path.λratio, λratio)
    push!(path.x, result.x)
    push!(path.converged, result.termination_status == MOI.OPTIMAL)
    push!(path.solve_time, result.solve_time)
    push!(path.node_count, result.node_count)
    push!(path.objective_value, result.objective_value)
    push!(path.datafit_value, value(f, A * result.x))
    push!(path.penalty_value, value(h, result.x))
    push!(path.support_size, norm(result.x, 0))
    push!(path.cv_mean, cv_mean)
    push!(path.cv_std, cv_std)
    return nothing
end

function compute_cv_statistics(x::Vector, f::AbstractDatafit, A::Matrix, nb_folds::Int)
    m, n = size(A)
    @assert nb_folds <= m
    m_fold = Int(ceil(m / nb_folds))
    cv_errors = Vector()
    for i = 1:nb_folds
        idx = randperm(m)[1:m_fold]
        f_idx = typeof(f)(f.y[idx])
        cv_error = value(f_idx, A[idx, :] * x)
        push!(cv_errors, cv_error)
    end
    return mean(cv_errors), std(cv_errors)
end

function isterminated(path::Path, options::PathOptions)
    (options.stop_if_unsolved & !path.converged[end]) && return true
    (path.support_size[end] > options.max_support_size) && return true
    return false
end

"""
    fit_path(
        solver::AbstractSolver,
        f::AbstractDatafit,
        h::AbstractPenalty,
        A::Matrix,
        kwargs...
    )

Fit a regularization [`Path`](@ref) for a [`Problem`](@ref). Additional keyword
arguments in `kwargs` are passed to a [`PathOptions`](@ref) instance.
"""
function fit_path(
    solver::AbstractSolver,
    f::AbstractDatafit,
    h::AbstractPenalty,
    A::Matrix;
    kwargs...,
)

    options = PathOptions(; kwargs...)

    @assert 0.0 < options.λratio_min <= options.λratio_max <= 1.0
    @assert 0 < options.λratio_num
    @assert options.nb_folds > 0

    λratio_sep =
        (log10(options.λratio_min) - log10(options.λratio_max)) / (options.λratio_num - 1)
    λratio_val = 10 .^ (log10(options.λratio_max):λratio_sep:log10(options.λratio_min))
    λmax = compute_λmax(f, h, A)
    x0 = zeros(size(A)[2])
    path = Path()

    options.verbosity && display_path_head()
    for λratio in λratio_val
        λ = λratio * λmax
        problem = Problem(f, h, A, λ)
        result = optimize(solver, problem, x0 = x0)
        fill_path!(path, f, h, A, λ, λratio, result, options)
        copy!(x0, result.x)
        options.verbosity && display_path_info(path)
        isterminated(path, options) && break
    end
    options.verbosity && display_path_tail()

    return path
end

function Base.show(io::IO, path::Path)
    display_path_head()
    for i = 1:length(path.λ)
        display_path_info(path, i)
    end
    display_path_tail()
end
