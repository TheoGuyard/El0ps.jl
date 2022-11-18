"""
    AbstractSolver

Abstract type for solver tailored to the problem min F(y,Ax) + λ * (||x||_0 + G(x)). 
"""
abstract type AbstractSolver end

"""
    AbstractResult

Abstract type to store results of an `AbstractSolver`. 
"""
abstract type AbstractResult end

"""
    optimize(
        solver::AbstractSolver,
        problem::Problem;
        x0::Union{Vector,Nothing}=nothing,
    )

Optimize `Problem` with the solver `solver` using the warm start `x0`. 
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