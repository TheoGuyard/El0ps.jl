abstract type AbstractDatafit end

function Base.show(io::IO, F::AbstractDatafit)
    error("Function 'show' is not implemented for the datafit $F")
end

function lipschitz_constant(F::AbstractDatafit, y::Vector)
    error("Function 'lipschitz_constant' is not implemented for the datafit $F")
end

function value(F::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'value' is not implemented for the datafit $F")
end

function gradient(F::AbstractDatafit, y::Vector, w::Vector)
    error("Function 'gradient' is not implemented for the datafit $F")
end

function conjugate(F::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'conjugate' is not implemented for the datafit $F")
end

function dual_scale!(F::AbstractDatafit, y::Vector, u::Vector)
    error("Function 'dual_scaling_factor' is not implemented for the datafit $F")
end
