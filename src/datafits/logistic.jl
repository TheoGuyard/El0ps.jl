"""
    Logistic <: AbstractDatafit

Logistic function \$f(\\mathbf{x}) = \\tfrac{1}{m} \\sum(\\log(\\mathbf{1} + \\exp(- \\mathbf{y} \\odot \\mathbf{x})))\$ where `m = length(y)`, where ``\\odot`` denotes the Hadamard product and where the `log` and the `exp` functions are taken component-wisely.
"""
struct Logistic <: AbstractDatafit 
    y::Vector
    m::Int
    Logistic(y::Vector) = new(y, length(y))
end

Base.show(io::IO, f::Logistic) = print(io, "Logistic")
dim_input(f::Logistic) = f.m
lipschitz_constant(f::Logistic) = 0.25 / f.m
value(f::Logistic, x::Vector) = sum(log.(1.0 .+ exp.(-f.y .* x))) / f.m
gradient(f::Logistic, x::Vector) = @. -f.y / (1. + exp(f.y * x)) / f.m
function conjugate(f::Logistic, x::Vector)
    v = -(x .* f.y) * f.m
    r = 1. .- v
    return (v' * log.(v) + r' * log.(r)) / f.m
end
