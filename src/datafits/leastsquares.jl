"""
    LeastSquares <: AbstractDatafit

Least-squares function.
"""
struct LeastSquares <: AbstractDatafit 
    y::Vector
    m::Int
end

"""
    LeastSquares(y::Vector)

Least-squares function `f(w) = norm(f.y - w, 2)^2 / f.m`.
"""
LeastSquares(y::Vector) = new(y, length(y))

Base.show(io::IO, f::LeastSquares) = print(io, "Least-Squares")
dim_input(f::LeastSquares) = f.m
lipschitz_constant(f::LeastSquares) = 1.0 / f.m
value(f::LeastSquares, x::Vector) = norm(x - f.y, 2)^2 / (2. * f.m)
gradient(f::LeastSquares, x::Vector) = (x - f.y) / f.m
conjugate(f::LeastSquares, x::Vector) = (0.5 * f.m) * (x' * x) + f.y' * x
