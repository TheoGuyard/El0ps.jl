"""
    AbstractPerturbation

Abstract type for the function `h` in a [`Problem`](@ref). This function is 
supposed to split. All its associated methods are therefore implemented with
one-dimensional inputs.
"""
abstract type AbstractPerturbation end

function Base.show(io::IO, h::AbstractPerturbation)
    error("Function 'show' is not implemented for the perturbation $h")
end

"""
    value_1d(h::AbstractPerturbation, x::Float64)

Value of `h`.
"""
function value_1d(h::AbstractPerturbation, x::Float64)
    error("Function 'value' is not implemented for the perturbation $h")
end

"""
    conjugate_1d(h::AbstractPerturbation, v::Float64)

Value of the conjugate of `h`.
"""
function conjugate_1d(h::AbstractPerturbation, v::Float64)
    error("Function 'conjugate_1d' is not implemented for the perturbation $h")
end

"""
    prox_1d(h::AbstractPerturbation, x::Float64, η::Float64)

Proximal mapping of the function `ηh`.
"""
function prox_1d(h::AbstractPerturbation, x::Float64, η::Float64)
    error("Function 'prox_1d' not implemented for the perturbation $h")
end

"""
    dual_scale!(h::AbstractPerturbation, A::Matrix, u::Vector, λ::Float64)

In-place transforms some vector `u` into a feasible one for the conjugate of `h` 
when evaluated at `A'u`. Returns the value of `A'u`.
"""
function dual_scale!(h::AbstractPerturbation, A::Matrix, u::Vector, λ::Float64)
    error("Function 'dual_scale!' is not implemented for the perturbation $h")
end

"""
    value(h::AbstractPerturbation, v::Vector)

Value of `h` evaluated with a vector `v` as input.
"""
function value(h::AbstractPerturbation, v::Vector)
    return length(v) > 0 ? sum(value_1d(h, vi) for vi in v) : 0.
end

"""
    conjugate_vectorized(h::AbstractPerturbation, v::Vector)

Coordinate-wise value of the conjugate of `h` when evaluated with a vector `v` 
as input.
"""
function conjugate_vectorized(h::AbstractPerturbation, v::Vector)
    return Vector{Float64}([conjugate_1d(h, vi) for vi in v])
end

"""
    conjugate(h::AbstractPerturbation, v::Vector)

Value of the conjugate of `h` evaluated with a vector `v` as input.
"""
function conjugate(h::AbstractPerturbation, v::Vector)
    return length(v) > 0 ? sum(conjugate_vectorized(h, v)) : 0.
end

"""
    prox(h::AbstractPerturbation, v::Vector, η::Float64)

Proximal mapping of the function `ηh` evaluated with a vector `v` as input.
"""
function prox(h::AbstractPerturbation, v::Vector, η::Float64)
    return Vector{Float64}([prox_1d(h, vi, η) for vi in v])
end

"""
    params_to_dict(h::AbstractPerturbation)

Returns a dictionary with the parameters name and value of the function `h`.
"""
function params_to_dict(h::AbstractPerturbation)
    htype = typeof(h)
    hfields = fieldnames(htype)
    params = OrderedDict(hfields .=> getfield.(Ref(h), hfields))
    filter!(kv -> !(first(kv) in [:μ,:τ]), params)
    return params
end

"""
    bind_model!(h::AbstractPerturbation, model::JuMP.Model)

Formulate the function `h` in `model`. See [`initialize_model`](@ref) for more 
details about the argument `model`.
"""
function bind_model!(h::AbstractPerturbation, model::JuMP.Model)
    error("Function 'bind_model!' is not implemented for the perturbation $h")
end
