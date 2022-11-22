"""
    AbstractDatafit

Abstract type for the function `F` in a [`Problem`](@ref).
"""
abstract type AbstractDatafit end

function Base.show(io::IO, F::AbstractDatafit)
    error("Function 'show' is not implemented for the datafit $F")
end

"""
    lipschitz_constant(F::AbstractDatafit, y::Vector)

Lischitz constant of the gradient of `F` with respect to its second argument.
"""
function lipschitz_constant(F::AbstractDatafit, y::Vector)
    error("Function 'lipschitz_constant' is not implemented for the datafit $F")
end

"""
    value(F::AbstractDatafit, y::Vector, w::Vector)

Value of `F`.
"""
function value(F::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'value' is not implemented for the datafit $F")
end

"""
    gradient(F::AbstractDatafit, y::Vector, w::Vector)

Gradient of `F` with respect to its second argument.
"""
function gradient(F::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'gradient' is not implemented for the datafit $F")
end

"""
    conjugate(F::AbstractDatafit, y::Vector, w::Vector)

Fenchel conjugate of `F` with respect to its second argument.
"""
function conjugate(F::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'conjugate' is not implemented for the datafit $F")
end

"""
    dual_scale!(F::AbstractDatafit, y::Vector, u::Vector)

In-place transforms a vector `u` into a feasible one for the conjugate of `F`.
"""
function dual_scale!(F::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'dual_scaling_factor' is not implemented for the datafit $F")
end

"""
    bind_model!(F::AbstractDatafit, y::Vector, model::JuMP.Model)

Formulate the function `F` in `model`. See [`initialize_model`](@ref) for more 
details about the argument `model`.
"""
function bind_model!(F::AbstractDatafit, y::Vector, model::JuMP.Model)
    error("Function 'bind_model!' is not implemented for the datafit $F")
end
