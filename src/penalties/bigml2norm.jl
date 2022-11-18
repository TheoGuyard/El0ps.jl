"""
    BigmL2norm

Big-M plus L2-norm function G(x) = (α/2) * ||x||_2^2 + Ind(||x||_Inf <= M).

# Arguments

- `M::Float64` : Big-M value.
- `α::Float64` : L2 regularization strength.
"""
struct BigmL2norm <: AbstractPenalty
    M::Float64
    α::Float64
    τ::Float64
    μ::Float64
    function BigmL2norm(M::Float64, α::Float64)
        (M > 0.) || error("Parameter M must be positive")
        (α > 0.) || error("Parameter α must be positive")
        τ = sqrt(2. / α) < M ? sqrt(2. * α) : (1. / M) + (α / 2.) * M
        μ = sqrt(2. / α) < M ? τ / α : M
        return new(M, α, τ, μ)
    end
end

Base.show(io::IO, G::BigmL2norm) = print(io, "Bigm + L2-norm")
value_1d(G::BigmL2norm, x::Float64) = abs(x) <= G.M ? (G.α / 2.) * x^2 : Inf
function conjugate_1d(G::BigmL2norm, v::Float64)
    r = clamp(v / G.α, -G.M, G.M)
    return v * r - (G.α / 2.) * r^2
end
prox_1d(G::BigmL2norm, x::Float64, η::Float64) = clamp(x / (1. + η * G.α), -G.M, G.M)
dual_scale!(G::BigmL2norm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(G::BigmL2norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, s[1:n] >= 0.0)
    for i in eachindex(s, model[:z], model[:x])
        @constraint(
            model,
            [0.5 * s[i]; model[:z][i]; model[:x][i]] in RotatedSecondOrderCone()
        )
    end
    @constraint(model, model[:x] .>= -G.M .* model[:z])
    @constraint(model, model[:x] .<= G.M .* model[:z])
    @constraint(model, model[:Ωcost] >= sum(model[:z]) + (G.α / 2.) * sum(s))
    return nothing
end
