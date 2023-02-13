"""
    L1norm

L1-norm function h(x) = α * ||x||_1.

# Arguments

- `α::Float64` : L1 regularization strength.
"""
struct L1norm <: AbstractPerturbation
    α::Float64
    τ::Float64
    μ::Float64
    function L1norm(α::Float64)
        (α > 0.) || error("Parameter α must be positive")
        return new(α, α, Inf)
    end
end

Base.show(io::IO, h::L1norm) = print(io, "L1-norm")
value_1d(h::L1norm, x::Float64) = h.α * abs(x)
conjugate_1d(h::L1norm, v::Float64) = (abs.(v) <= h.α) ? 0. : Inf
prox_1d(h::L1norm, x::Float64, η::Float64) = sign(x) * max(abs(x) - η * h.α, 0.0)

function dual_scale!(h::L1norm, A::Matrix, u::Vector, λ::Float64)
    v = A' * u
    s = λ * h.α / norm(v, Inf)
    lmul!(s, u)
    lmul!(s, v)
    return v
end

function bind_model!(h::L1norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, xabs[1:n])
    @constraint(model, xabs .* model[:z] .>= model[:x])
    @constraint(model, xabs .* model[:z] .>= -model[:x])
    @constraint(model, model[:Ωcost] >= sum(model[:z]) + h.α * sum(xabs))
    return nothing
end
