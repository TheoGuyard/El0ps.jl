struct LeastSquares <: AbstractDatafit end

Base.show(io::IO, F::LeastSquares) = print(io, "Least-Squares")
lipschitz_constant(F::LeastSquares, y::Vector) = 1.0 / length(y)

function value(F::LeastSquares, y::Vector, w::Vector)
    return norm(w - y, 2)^2 / (2. * length(y))
end

function gradient(F::LeastSquares, y::Vector, w::Vector)
    return (w - y) ./ length(y)
end

function conjugate(F::LeastSquares, y::Vector, u::Vector)
    return (0.5 * length(y)) * norm(u, 2)^2 + y' * u
end

function dual_scale!(F::LeastSquares, y::Vector, u::Vector)
    return nothing
end
