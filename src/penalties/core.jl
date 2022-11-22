"""
    AbstractPenalty

Abstract type for the function `G` in a [`Problem`](@ref). This function is 
supposed to split. All its associated methods are therefore implemented with
one-dimensional inputs.
"""
abstract type AbstractPenalty end

function Base.show(io::IO, G::AbstractPenalty)
    error("Function 'show' is not implemented for the penalty $G")
end

"""
    value_1d(G::AbstractPenalty, x::Float64)

Value of `G`.
"""
function value_1d(G::AbstractPenalty, x::Float64)
    error("Function 'value' is not implemented for the penalty $G")
end

"""
    conjugate_1d(G::AbstractPenalty, v::Float64)

Value of the conjugate of `G`.
"""
function conjugate_1d(G::AbstractPenalty, v::Float64)
    error("Function 'conjugate_1d' is not implemented for the penalty $G")
end

"""
    prox_1d(G::AbstractPenalty, x::Float64, η::Float64)

Proximal mapping of the function `ηG`.
"""
function prox_1d(G::AbstractPenalty, x::Float64, η::Float64)
    error("Function 'prox_1d' not implemented for the penalty $G")
end

"""
    dual_scale!(G::AbstractPenalty, A::Matrix, u::Vector, λ::Float64)

In-place transforms some vector `u` into a feasible one for the conjugate of `G` 
when evaluated at `A'u`. Returns the value of `A'u`.
"""
function dual_scale!(G::AbstractPenalty, A::Matrix, u::Vector, λ::Float64)
    error("Function 'dual_scale!' is not implemented for the penalty $G")
end

"""
    value(G::AbstractPenalty, v::Vector)

Value of `G` evaluated with a vector `v` as input.
"""
function value(G::AbstractPenalty, v::Vector)
    return length(v) > 0 ? sum(value_1d(G, vi) for vi in v) : 0.
end

"""
    conjugate_vectorized(G::AbstractPenalty, v::Vector)

Coordinate-wise value of the conjugate of `G` when evaluated with a vector `v` 
as input.
"""
function conjugate_vectorized(G::AbstractPenalty, v::Vector)
    return Vector{Float64}([conjugate_1d(G, vi) for vi in v])
end

"""
    conjugate(G::AbstractPenalty, v::Vector)

Value of the conjugate of `G` evaluated with a vector `v` as input.
"""
function conjugate(G::AbstractPenalty, v::Vector)
    return length(v) > 0 ? sum(conjugate_vectorized(G, v)) : 0.
end

"""
    prox(G::AbstractPenalty, v::Vector, η::Float64)

Proximal mapping of the function `ηG` evaluated with a vector `v` as input.
"""
function prox(G::AbstractPenalty, v::Vector, η::Float64)
    return Vector{Float64}([prox_1d(G, vi, η) for vi in v])
end

"""
    params_to_dict(G::AbstractPenalty)

Returns a dictionary with the parameters names of the function G(.) and their 
value. The dictionary is empty is G(.) has no parameters.
"""
function params_to_dict(G::AbstractPenalty)
    Gtype = typeof(G)
    Gfields = fieldnames(Gtype)
    params = Dict(Gfields .=> getfield.(Ref(G), Gfields))
    filter!(kv -> !(first(kv) in [:μ,:τ]), params)
    return params
end

"""
    bind_model!(G::AbstractPenalty, model::JuMP.Model)

Formulate the function G(.) in `model` using an epigraph method. See
`initialize_model` for more details.
"""
function bind_model!(G::AbstractPenalty, model::JuMP.Model)
    error("Function 'bind_model!' is not implemented for the penalty $G")
end
