"""
    BoundingType

Type of the bounding operation.
"""
@enum BoundingType begin
    LOWER_BOUNDING
    UPPER_BOUNDING
end

"""
    AbstractBoundingSolver

Abstract supertype for the boudning solver in a [`BnbSolver`](@ref).
"""
abstract type AbstractBoundingSolver end

"""
    bounding_type(bounding_solver::AbstractBoundingSolver)

Bounding type of the `bounding_solver`.
"""
bounding_type(bounding_solver::AbstractBoundingSolver) = error("Not implemented")

"""
    bound!(
        bounding_solver::AbstractBoundingSolver, 
        problem::Problem, 
        solver,
        node,
        options,
    )

Perform the bounding operation at a given `node` of a Branch-and-Bound `solver`
during the solving process of a [`Problem`](@ref).
"""
bound!(bounding_solver::AbstractBoundingSolver, problem::Problem, solver, node, options) =
    error("Not implemented")
