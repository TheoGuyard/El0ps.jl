"""
    AbstractSolver

Abstract type for a solver tailored to a [`Problem`](@ref).
"""
abstract type AbstractSolver end

"""
    AbstractResult

Abstract type to store results of an [`AbstractSolver`](@ref). 
"""
abstract type AbstractResult end

"""
    optimize(
        solver::AbstractSolver,
        problem::Problem;
        x0::Union{Vector,Nothing}=nothing,
    )

Optimize a [`Problem`](@ref) with an [`AbstractSolver`](@ref). The argument `x0` 
is used as a warm start.  
"""
function optimize(
    solver::AbstractSolver,
    problem::Problem;
    x0::Union{Vector,Nothing}=nothing,
    )
    error("Function 'optimize' not implemented for solver $solver")
end

function Base.print(io::IO, result::AbstractResult)
    println(io, "Result")
    println(io, "  Status     : $(result.termination_status)")
    println(io, "  Objective  : $(result.objective_value)")
    println(io, "  Non-zeros  : $(norm(result.x, 0))")
    println(io, "  Last gap   : $(result.relative_gap)")
    println(io, "  Solve time : $(result.solve_time) seconds")
    print(io, "  Node count : $(result.node_count)")
end
