"""
    Logistic <: AbstractDatafit

Logistic function `f(w) = sum(log.(1.0 .+ exp.(-f.y .* x))) / f.m` where 
`f.m = length(f.y)`.

# Attributes 

- `y::Vector` : Target data vector.
- `m::Int` : Size of `y`.
"""
struct Logistic <: AbstractDatafit
    y::Vector
    m::Int
end


"""
    Logistic(y::Vector)

[`Logistic`](@ref) constructor.
"""
Logistic(y::Vector) = Logistic(y, length(y))

Base.show(io::IO, f::Logistic) = print(io, "Logistic")
dim_input(f::Logistic) = f.m
lipschitz_constant(f::Logistic) = 0.25 / f.m
value(f::Logistic, x::Vector) = sum(log.(1.0 .+ exp.(-f.y .* x))) / f.m
gradient(f::Logistic, x::Vector) = @. -f.y / (1.0 + exp(f.y * x)) / f.m
function conjugate(f::Logistic, x::Vector)
    v = -(x .* f.y) * f.m
    r = 1.0 .- v
    return (v' * log.(v) + r' * log.(r)) / f.m
end
