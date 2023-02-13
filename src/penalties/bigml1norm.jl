"""
    BigmL1norm

Big-M plus L1-norm function h(x) = α||x||_1 + Ind(||x||_Inf <= M).

# Arguments

- `M::Float64` : Big-M value.
- `α::Float64` : L1 regularization strength.
"""
struct BigmL1norm <: AbstractPerturbation
    M::Float64
    α::Float64
    τ::Float64
    μ::Float64
    function BigmL1norm(M::Float64, α::Float64)
        (M > 0.) || error("Parameter M must be positive")
        (α > 0.) || error("Parameter α must be positive")
        return new(M, α, (1. / M) + α, M)
    end
end

Base.show(io::IO, h::BigmL1norm) = print(io, "Bigm + L1-norm")
value_1d(h::BigmL1norm, x::Float64) = abs(x) <= h.M ? h.α * abs(x) : Inf
conjugate_1d(h::BigmL1norm, v::Float64) = h.M * (max(abs(v) - h.α, 0.))
prox_1d(h::BigmL1norm, x::Float64, η::Float64) = sign(x) * clamp(abs(x) - η * h.α, 0., h.M)
dual_scale!(h::BigmL1norm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(h::BigmL1norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, xabs[1:n])
    @constraint(model, xabs .>= model[:x])
    @constraint(model, xabs .>= -model[:x])
    @constraint(model, model[:x] .>= -h.M .* model[:z])
    @constraint(model, model[:x] .<= h.M .* model[:z])
    @constraint(model, model[:Gcost] >= sum(model[:z]) + h.α * sum(model[:xabs]))
    return nothing
end
