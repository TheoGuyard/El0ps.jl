"""
    LeastSquares <: AbstractDatafit

Least-squares function 

``f(\\mathbf{y},\\mathbf{w}) = \\tfrac{1}{m} \\|\\mathbf{y}-\\mathbf{w}\\|_2^2``

where `m = length(y)`.
"""
struct LeastSquares <: AbstractDatafit end

Base.show(io::IO, f::LeastSquares) = print(io, "Least-Squares")

function lipschitz_constant(f::LeastSquares, y::Vector)
    return 1.0 / length(y)
end

function value(f::LeastSquares, y::Vector, w::Vector)
    return norm(w - y, 2)^2 / (2. * length(y))
end

function gradient(f::LeastSquares, y::Vector, w::Vector)
    return (w - y) ./ length(y)
end

function conjugate(f::LeastSquares, y::Vector, u::Vector)
    return (0.5 * length(y)) * norm(u, 2)^2 + y' * u
end

function dual_scale!(f::LeastSquares, y::Vector, u::Vector)
    return nothing
end

function bind_model!(f::LeastSquares, y::Vector, model::JuMP.Model)
    @constraint(
        model, 
        model[:fcost] >= 0.5 * (y - model[:w])' * (y - model[:w]) / length(y)
    )
    return nothing
end
