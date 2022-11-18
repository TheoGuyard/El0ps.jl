"""
    Logistic

Logistic function F(y,w) = (1/m) * sum(log(1 + exp(- y ⊙ w))) where m=size(y),
where ⊙ denotes the Hadamard product and where the `log` and the `exp` 
functions are taken component-wisely.
"""
struct Logistic <: AbstractDatafit end

Base.show(io::IO, F::Logistic) = print(io, "Logistic")
lipschitz_constant(F::Logistic, y::Vector) = 0.25 / length(y)

function value(F::Logistic, y::Vector, w::Vector)
    return sum(log.(1.0 .+ exp.(-y .* w))) / length(y)
end

function gradient(F::Logistic, y::Vector, w::Vector)
    return (@. -y / (1. + exp(y * w))) / length(y)
end

function conjugate(F::Logistic, y::Vector, u::Vector)
    v = -(u .* y) * length(y)
    r = 1. .- v
    return (v' * log.(v) + r' * log.(r)) / length(y)
end

function dual_scale!(F::Logistic, y::Vector, u::Vector)
    m = length(y)
    for j in eachindex(y, u)
        u[j] = clamp(u[j], min(y[j] / m, 0.), max(y[j] / m, 0.))
    end
    return nothing
end

function bind_model!(F::Logistic, y::Vector, model::JuMP.Model)
    m = length(y)
    @variable(model, l[1:m])
    @variable(model, u[1:m] >= 0)
    @variable(model, v[1:m] >= 0)
    for j in eachindex(y, l, u, v, model[:w])
        @constraint(model, 1. >= u[j] + v[j])
        @constraint(model, [y[j] * model[:w][j] - l[j]; 1.; u[j]] in MOI.ExponentialCone())
        @constraint(model, [-l[j]; 1.; v[j]] in MOI.ExponentialCone())
    end
    @constraint(model, model[:Fcost] >= sum(l) / m)
    return nothing
end
