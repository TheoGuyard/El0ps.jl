abstract type AbstractBoundingSolver end

@enum BoundingType begin
    LOWER
    UPPER
end

function bound!(
    bounding_solver::AbstractBoundingSolver, 
    problem, 
    solver, 
    node, 
    options,   
    bounding_type::BoundingType,
    )
    error("Function 'bound!' not implemented for the solver '$solver' and the bounding type '$bounding_type'")
end

function prox_l1_1d(x::Float64, η::Float64)
    return sign(x) * max(abs(x) - η, 0.0)
end

function prox_l1(x::Vector{Float64}, η::Float64)
    return @. sign(x) * max(abs(x) - η, 0.0)
end