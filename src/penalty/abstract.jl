"""
    AbstractPenalty

Abstract supertype for the function `h` in a [`Problem`](@ref).
"""
abstract type AbstractPenalty end

"""
    compute_τ(h::AbstractPenalty, λ::Float64)

Compute the `τ` value of the funtion `g(x) = λ * norm(x, 0) + h(x)`.
"""
compute_τ(h::AbstractPenalty, λ::Float64) = error("Not implemented")

"""
    compute_μ(h::AbstractPenalty, λ::Float64)

Compute the `μ` value of the funtion `g(x) = λ * norm(x, 0) + h(x)`.
"""
compute_μ(h::AbstractPenalty, λ::Float64) = error("Not implemented")


"""
    value_1d(h::AbstractPenalty, x::Float64)

Value of a splitting term of `h` at `x`.
"""
value_1d(h::AbstractPenalty, x::Float64) = error("Not implemented")

"""
    value(h::AbstractPenalty, x::Vector)

Value of `h` at `x`.
"""
value(h::AbstractPenalty, x::Vector) = sum(value_1d(h, xi) for xi in x)


"""
    conjugate_1d(h::AbstractPenalty, x::Float64)

Value of the conjugate of a splitting term of `h` at `x`.
"""
conjugate_1d(h::AbstractPenalty, x::Float64) = error("Not implemented")

"""
    conjugate(h::AbstractPenalty, x::Vector)

Value of the conjugate of `h` at `x`.
"""
conjugate(h::AbstractPenalty, x::Vector) = sum(conjugate_1d(h, xi) for xi in x)


"""
    prox_1d(h::AbstractPenalty, x::Float64, η::Float64)

Proximal mapping of a splitting term of `ηh` at `x`.
"""
prox_1d(h::AbstractPenalty, x::Float64, η::Float64) = error("Not implemented")

"""
    prox(h::AbstractPenalty, x::Vector, η::Float64)

Proximal mapping of `ηh` evaluated at `x`.
"""
prox(h::AbstractPenalty, x::Vector, η::Float64) = [prox_1d(h, xi, η) for xi in x]

"""
    params_to_dict(h::AbstractPenalty)

Returns a dictionary with the parameters name and value of the function `h`.
"""
function params_to_dict(h::AbstractPenalty)
    htype = typeof(h)
    hfields = fieldnames(htype)
    params = OrderedDict(hfields .=> getfield.(Ref(h), hfields))
    filter!(kv -> !(first(kv) in [:μ, :τ]), params)
    return params
end

"""
    prox_l1_1d(x::Float64, η::Float64)

Proximal mapping of the function `η|.|` at `x`.
"""
prox_l1_1d(x::Float64, η::Float64) = sign(x) * max(abs(x) - η, 0.0)
