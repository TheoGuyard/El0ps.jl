"""
    BigmL2norm

Big-M plus L2-norm function h(x) = α * ||x||_2^2 + Ind(||x||_Inf <= M).

# Arguments

- `M::Float64` : Big-M value.
- `α::Float64` : L2 regularization strength.
"""
struct BigmL2norm <: AbstractPerturbation
    M::Float64
    α::Float64
    τ::Float64
    μ::Float64
    function BigmL2norm(M::Float64, α::Float64)
        (M > 0.) || error("Parameter M must be positive")
        (α > 0.) || error("Parameter α must be positive")
        τ = sqrt(1. / α) < M ? sqrt(4. * α) : (1. / M) + α * M
        μ = sqrt(1. / α) < M ? τ / (2. * α) : M
        return new(M, α, τ, μ)
    end
end

Base.show(io::IO, h::BigmL2norm) = print(io, "Bigm + L2-norm")
value_1d(h::BigmL2norm, x::Float64) = abs(x) <= h.M ? h.α * x^2 : Inf
function conjugate_1d(h::BigmL2norm, v::Float64)
    r = clamp(v / (2. * h.α), -h.M, h.M)
    return v * r - h.α * r^2
end
prox_1d(h::BigmL2norm, x::Float64, η::Float64) = clamp(x / (1. + 2. * η * h.α), -h.M, h.M)
dual_scale!(h::BigmL2norm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(h::BigmL2norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, s[1:n] >= 0.0)
    for i in eachindex(s, model[:z], model[:x])
        @constraint(
            model,
            [0.5 * s[i]; model[:z][i]; model[:x][i]] in RotatedSecondOrderCone()
        )
    end
    @constraint(model, model[:x] .>= -h.M .* model[:z])
    @constraint(model, model[:x] .<= h.M .* model[:z])
    @constraint(model, model[:Gcost] >= sum(model[:z]) + h.α * sum(s))
    return nothing
end
