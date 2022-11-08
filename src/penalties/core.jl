abstract type AbstractPenalty end

function Base.show(io::IO, G::AbstractPenalty)
    error("Function 'show' is not implemented for the penalty $G")
end

function value_1d(G::AbstractPenalty, x::Float64)
    error("Function 'value' is not implemented for the penalty $G")
end

function conjugate_1d(G::AbstractPenalty, v::Float64)
    error("Function 'conjugate_1d' is not implemented for the penalty $G")
end

function prox_1d(G::AbstractPenalty, x::Float64, η::Float64)
    error("Function 'prox_1d' not implemented for the penalty $G")
end

function dual_scale!(G::AbstractPenalty, A::Matrix, u::Vector, λ::Float64)
    error("Function 'dual_scale!' is not implemented for the penalty $G")
end

# TODO : optimize with vectorization ?
value(G::AbstractPenalty, v::Vector) = length(v) > 0 ? sum(value_1d(G, vi) for vi in v) : 0.
conjugate(G::AbstractPenalty, v::Vector) = length(v) > 0 ? sum(conjugate_1d(G, vi) for vi in v) : 0.
conjugate_vectorized(G::AbstractPenalty, v::Vector) = [conjugate_1d(G, vi) for vi in v]
prox(G::AbstractPenalty, x::Vector, η::Float64) = [prox_1d(G, xi, η) for xi in x]

function params_to_dict(G::AbstractPenalty)
    Gtype = typeof(G)
    Gfields = fieldnames(Gtype)
    params = Dict(Gfields .=> getfield.(Ref(G), Gfields))
    filter!(kv -> !(first(kv) in [:μ,:τ]), params)
    return params
end

function bind_model!(G::AbstractPenalty, model::JuMP.Model)
    error("Function 'bind_model!' is not implemented for the penalty $G")
end

