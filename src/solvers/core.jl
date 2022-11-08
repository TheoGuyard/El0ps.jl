abstract type AbstractSolver end
abstract type AbstractOptions end
abstract type AbstractResult end

function optimize(
    solver::AbstractSolver,
    problem::Problem;
    x0::Union{Vector{Float64},Nothing}=nothing,
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