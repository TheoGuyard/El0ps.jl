"""
    AbstractPenalty

Abstract type for the function G(.) in `Problem`.
"""
abstract type AbstractPenalty end

function Base.show(io::IO, G::AbstractPenalty)
    error("Function 'show' is not implemented for the penalty $G")
end

"""
    value_1d(G::AbstractPenalty, x::Float64)

Value of G(.) evaluated at some scalar x.
"""
function value_1d(G::AbstractPenalty, x::Float64)
    error("Function 'value' is not implemented for the penalty $G")
end

"""
    conjugate_1d(G::AbstractPenalty, v::Float64)

Conjugate of the function G(.) evaluated at some scalar v.
"""
function conjugate_1d(G::AbstractPenalty, v::Float64)
    error("Function 'conjugate_1d' is not implemented for the penalty $G")
end

"""
    prox_1d(G::AbstractPenalty, x::Float64, η::Float64)

Proximal mapping of the function ηG(.) evaluated at some scalar x.
"""
function prox_1d(G::AbstractPenalty, x::Float64, η::Float64)
    error("Function 'prox_1d' not implemented for the penalty $G")
end

"""
    dual_scale!(G::AbstractPenalty, A::Matrix, u::Vector, λ::Float64)

In-place transform any vector u into a feasible one for the conjugate of the 
function G(.) evaluated at A'u. Returns the value of A'u.
"""
function dual_scale!(G::AbstractPenalty, A::Matrix, u::Vector, λ::Float64)
    error("Function 'dual_scale!' is not implemented for the penalty $G")
end

# TODO : optimize with vectorization ?
value(G::AbstractPenalty, v::Vector) = length(v) > 0 ? sum(value_1d(G, vi) for vi in v) : 0.
conjugate(G::AbstractPenalty, v::Vector) = length(v) > 0 ? sum(conjugate_1d(G, vi) for vi in v) : 0.
conjugate_vectorized(G::AbstractPenalty, v::Vector) = [conjugate_1d(G, vi) for vi in v]
prox(G::AbstractPenalty, x::Vector, η::Float64) = [prox_1d(G, xi, η) for xi in x]

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

