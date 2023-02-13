"""
    L2norm

L2-norm function h(x) = α * ||x||_2^2.

# Arguments

- `α::Float64` : L2 regularization strength.
"""
struct L2norm <: AbstractPerturbation
    α::Float64
    τ::Float64
    μ::Float64
    function L2norm(α::Float64)
        (α > 0.) || error("Parameter α must be positive")
        τ = 2. * sqrt(α)
        μ = 1. / sqrt(α)
        return new(α, τ, μ)
    end
end

Base.show(io::IO, h::L2norm) = print(io, "L2-norm")
value_1d(h::L2norm, x::Float64) = h.α * x^2
conjugate_1d(h::L2norm, v::Float64) = v^2 / (4. * h.α)
prox_1d(h::L2norm, x::Float64, η::Float64) = x / (1. + 2. * η * h.α)
dual_scale!(h::L2norm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(h::L2norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, s[1:n] >= 0.0)
    for i in eachindex(s, model[:z], model[:x])
        @constraint(
            model,
            [0.5 * s[i]; model[:z][i]; model[:x][i]] in RotatedSecondOrderCone()
        )
    end
    @constraint(model, model[:Ωcost] >= sum(model[:z]) + h.α *  sum(s))
    return nothing
end
