"""
    L1L2norm <: AbstractPenalty

L1L2-norm function `h(x) = h.α * norm(x, 1) + h.β * norm(x, 2)^2`, where
`h.α > 0` and `h.β > 0`.

# Arguments

- `α::Float64` : L1 regularization strength.
- `β::Float64` : L2 regularization strength.
"""
struct L1L2norm <: AbstractPenalty
    α::Float64
    β::Float64
    function L1L2norm(α::Float64, β::Float64)
        (α > 0.0) || error("Parameter α must be positive")
        (β > 0.0) || error("Parameter β must be positive")
        return new(α, β)
    end
end

Base.show(io::IO, h::L1L2norm) = print(io, "L1L2-norm")
compute_τ(h::L1L2norm, λ::Float64) = h.α + sqrt(4.0 * h.β * λ)
compute_μ(h::L1L2norm, λ::Float64) = sqrt(λ / h.β)
function compute_λmax(f::AbstractDatafit, h::L1L2norm, A::Matrix)
    v = norm(A' * gradient(f, zeros(dim_input(f))), Inf)
    return max((v - h.α)^2 / (4.0 * h.β), 0.0)
end
value_1d(h::L1L2norm, x::Float64) = h.α * abs(x) + h.β * x^2
conjugate_1d(h::L1L2norm, v::Float64) = max(abs(v) - h.α, 0.0)^2 / (4.0 * h.β)
function prox_1d(h::L1L2norm, x::Float64, η::Float64)
    return (sign(x) / (1.0 + 2.0 * η * h.β)) * max(abs(x) - η * h.α, 0.0)
end
