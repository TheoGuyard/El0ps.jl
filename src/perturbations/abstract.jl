"""
    AbstractPerturbation

Abstract supertype for the function `h` in a [`Problem`](@ref).
"""
abstract type AbstractPerturbation end

"""
    value_1d(h::AbstractPerturbation, x::Float64)

Value of a splitting term of `h` at `x`.
"""
value_1d(h::AbstractPerturbation, x::Float64) = error("Not implemented")

"""
    value(h::AbstractPerturbation, x::Vector)

Value of `h` at `x`.
"""
value(h::AbstractPerturbation, x::Vector) = sum(value_1d(h, xi) for xi in x)


"""
    conjugate_1d(h::AbstractPerturbation, x::Float64)

Value of the conjugate of a splitting term of `h` at `x`.
"""
conjugate_1d(h::AbstractPerturbation, x::Float64) = error("Not implemented")

"""
    conjugate(h::AbstractPerturbation, x::Vector)

Value of the conjugate of `h` at `x`.
"""
conjugate(h::AbstractPerturbation, x::Vector) = sum(conjugate_1d(h, xi) for xi in x)


"""
    prox_1d(h::AbstractPerturbation, x::Float64, η::Float64)

Proximal mapping of a splitting term of `ηh` at `x`.
"""
prox_1d(h::AbstractPerturbation, x::Float64, η::Float64) = error("Not implemented")

"""
    prox(h::AbstractPerturbation, x::Vector, η::Float64)

Proximal mapping of `ηh` evaluated at `x`.
"""
prox(h::AbstractPerturbation, x::Vector, η::Float64) = [prox_1d(h, xi, η) for xi in x]

"""
    params_to_dict(h::AbstractPerturbation)

Returns a dictionary with the parameters name and value of the function `h`.
"""
function params_to_dict(h::AbstractPerturbation)
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
