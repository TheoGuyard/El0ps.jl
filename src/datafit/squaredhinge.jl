"""
    SquaredHinge <: AbstractDatafit

Squared-Hinge function `f(w) = norm(max.(1.- y .* w, 0.), 2)^2 / f.m` where
`f.m = length(f.y)`.

# Attributes

- `y::Vector` : Target data vector.
- `m::Int` : Size of `y`.
"""
struct SquaredHinge <: AbstractDatafit
    y::Vector
    m::Int
end

"""
    SquaredHinge(y::Vector)

[`SquaredHinge`](@ref) constructor.
"""
function SquaredHinge(y::Vector)
    @assert all(yi in [-1.0, 1.0] for yi in y)
    return SquaredHinge(y, length(y))
end

Base.show(io::IO, f::SquaredHinge) = print(io, "Squared-Hinge")
dim_input(f::SquaredHinge) = f.m
lipschitz_constant(f::SquaredHinge) = 2.0 / f.m
value(f::SquaredHinge, x::Vector) = norm(max.(1.0 .- f.y .* x, 0.0), 2)^2 / f.m
gradient(f::SquaredHinge, x::Vector) = -f.y .* max.(1.0 .- f.y .* x, 0.0) / (0.5 * f.m)
conjugate(f::SquaredHinge, x::Vector) =
    (0.5 * f.m) * (x' * x) + f.y' * x - norm(max.(-0.5 * f.m * (f.y .* x), 0.0), 2)^2 / f.m
