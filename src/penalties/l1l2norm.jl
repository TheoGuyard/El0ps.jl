"""
    L1L2norm

L1L2-norm function G(x) = α * ||x||_1 + β * ||x||_2^2.

# Arguments

- `α::Float64` : L1 regularization strength.
- `β::Float64` : L2 regularization strength.
"""
struct L1L2norm <: AbstractPenalty
    α::Float64
    β::Float64
    τ::Float64
    μ::Float64
    function L1L2norm(α::Float64, β::Float64)
        (α > 0.) || error("Parameter α must be positive")
        (β > 0.) || error("Parameter β must be positive")
        τ =  α + sqrt(4. * β)
        μ = sqrt(1. / β)
        return new(α, β, τ, μ)
    end
end

Base.show(io::IO, G::L1L2norm) = print(io, "L1L2-norm")
value_1d(G::L1L2norm, x::Float64) = G.α * abs(x) + G.β * x^2
conjugate_1d(G::L1L2norm, v::Float64) = max(abs(v) - G.α, 0.)^2 / (4. * G.β)
prox_1d(G::L1L2norm, x::Float64, η::Float64) = (sign(x) / (1. + 2. * η * G.β)) * max(abs(x) - η * G.α, 0.)
dual_scale!(G::L1L2norm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(G::L1L2norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, xabs[1:n])
    @variable(model, s[1:n] >= 0.0)
    @constraint(model, xabs .>= model[:x])
    @constraint(model, xabs .>= -model[:x])
    for i in eachindex(s, model[:z], model[:x])
        @constraint(
            model,
            [0.5 * s[i]; model[:z][i]; model[:x][i]] in RotatedSecondOrderCone()
        )
    end
    @constraint(model, model[:Ωcost] >= sum(model[:z]) + G.α *  sum(xabs) + G.β *  sum(s))
    return nothing
end
