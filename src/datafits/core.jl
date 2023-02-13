"""
    AbstractDatafit

Abstract type for the function `f` in a [`Problem`](@ref).
"""
abstract type AbstractDatafit end

function Base.show(io::IO, f::AbstractDatafit)
    error("Function 'show' is not implemented for the datafit $f")
end

"""
    lipschitz_constant(f::AbstractDatafit, y::Vector)

Lischitz constant of the gradient of `f` with respect to its second argument.
"""
function lipschitz_constant(f::AbstractDatafit, y::Vector)
    error("Function 'lipschitz_constant' is not implemented for the datafit $f")
end

"""
    value(f::AbstractDatafit, y::Vector, w::Vector)

Value of `f`.
"""
function value(f::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'value' is not implemented for the datafit $f")
end

"""
    gradient(f::AbstractDatafit, y::Vector, w::Vector)

Gradient of `f` with respect to its second argument.
"""
function gradient(f::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'gradient' is not implemented for the datafit $f")
end

"""
    conjugate(f::AbstractDatafit, y::Vector, w::Vector)

Fenchel conjugate of `f` with respect to its second argument.
"""
function conjugate(f::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'conjugate' is not implemented for the datafit $f")
end

"""
    dual_scale!(f::AbstractDatafit, y::Vector, u::Vector)

In-place transforms a vector `u` into a feasible one for the conjugate of `f`.
"""
function dual_scale!(f::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'dual_scaling_factor' is not implemented for the datafit $f")
end

"""
    bind_model!(f::AbstractDatafit, y::Vector, model::JuMP.Model)

Formulate the function `f` in `model`. See [`initialize_model`](@ref) for more 
details about the argument `model`.
"""
function bind_model!(f::AbstractDatafit, y::Vector, model::JuMP.Model)
    error("Function 'bind_model!' is not implemented for the datafit $f")
end
