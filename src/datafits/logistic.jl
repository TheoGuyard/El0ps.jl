"""
    Logistic <: AbstractDatafit

Logistic function 
    
``f(\\mathbf{y},\\mathbf{w}) = \\tfrac{1}{m} \\sum(\\log(\\mathbf{1} + \\exp(- \\mathbf{y} \\odot \\mathbf{w})))``

where `m = length(y)`, where ``\\odot`` denotes the Hadamard product and where 
the `log` and the `exp` functions are taken component-wisely.
"""
struct Logistic <: AbstractDatafit end

Base.show(io::IO, f::Logistic) = print(io, "Logistic")

function lipschitz_constant(f::Logistic, y::Vector)
    return 0.25 / length(y)
end

function value(f::Logistic, y::Vector, w::Vector)
    return sum(log.(1.0 .+ exp.(-y .* w))) / length(y)
end

function gradient(f::Logistic, y::Vector, w::Vector)
    return (@. -y / (1. + exp(y * w))) / length(y)
end

function conjugate(f::Logistic, y::Vector, u::Vector)
    v = -(u .* y) * length(y)
    r = 1. .- v
    return (v' * log.(v) + r' * log.(r)) / length(y)
end

function dual_scale!(f::Logistic, y::Vector, u::Vector)
    m = length(y)
    for j in eachindex(y, u)
        u[j] = clamp(u[j], min(y[j] / m, 0.), max(y[j] / m, 0.))
    end
    return nothing
end

function bind_model!(f::Logistic, y::Vector, model::JuMP.Model)
    m = length(y)
    @variable(model, l[1:m])
    @variable(model, u[1:m] >= 0)
    @variable(model, v[1:m] >= 0)
    for j in eachindex(y, l, u, v, model[:w])
        @constraint(model, 1. >= u[j] + v[j])
        @constraint(
            model, 
            [y[j] * model[:w][j] - l[j]; 1.; u[j]] in MOI.ExponentialCone()
        )
        @constraint(model, [-l[j]; 1.; v[j]] in MOI.ExponentialCone())
    end
    @constraint(model, model[:fcost] >= sum(l) / m)
    return nothing
end
