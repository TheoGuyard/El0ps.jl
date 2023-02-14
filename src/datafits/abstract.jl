"""
    AbstractDatafit

Abstract supertype for the function `f` in a [`Problem`](@ref).
"""
abstract type AbstractDatafit end

"""
    dim_input(f::AbstractDatafit)

Dimension required for the vector input of `f`.
"""
dim_input(f::AbstractDatafit) = error("Not implemented")

"""
    lipschitz_constant(f::AbstractDatafit)

Lischitz constant of the gradient of `f`.
"""
lipschitz_constant(f::AbstractDatafit) = error("Not implemented")

"""
    value(f::AbstractDatafit, x::Vector)

Value of `f` at `x`.
"""
value(f::AbstractDatafit, x::Vector) = error("Not implemented")

"""
    gradient(f::AbstractDatafit, x::Vector)

Gradient of `f` at `x`.
"""
gradient(f::AbstractDatafit, x::Vector) = error("Not implemented")

"""
    conjugate(f::AbstractDatafit, x::Vector)

Value of the conjugate of `f` at `x`.
"""
conjugate(f::AbstractDatafit, x::Vector) = error("Not implemented")

"""
    params_to_dict(f::AbstractDatafit)

Returns a dictionary with the parameters name and value of `f`.
"""
function params_to_dict(f::AbstractDatafit)
    ftype = typeof(f)
    ffields = fieldnames(ftype)
    params = OrderedDict(ffields .=> getfield.(Ref(f), ffields))
    return params
end
