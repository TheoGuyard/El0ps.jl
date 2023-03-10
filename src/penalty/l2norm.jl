"""
    L2norm <: AbstractPenalty

L2-norm function `h(x) = h.α * norm(x, 2)^2`, where `h.α > 0`.

# Arguments

- `α::Float64` : L2 regularization strength.
"""
struct L2norm <: AbstractPenalty
    α::Float64
    function L2norm(α::Float64)
        (α > 0.0) || error("Parameter α must be positive")
        return new(α)
    end
end

Base.show(io::IO, h::L2norm) = print(io, "L2-norm")
compute_τ(h::L2norm, λ::Float64) = sqrt(4.0 * h.α * λ)
compute_μ(h::L2norm, λ::Float64) = sqrt(λ / h.α)
function compute_λmax(f::AbstractDatafit, h::L2norm, A::Matrix)
    v = norm(A' * gradient(f, zeros(dim_input(f))), Inf)
    return max(v^2 / (4. * h.α), 0.)
end
value_1d(h::L2norm, x::Float64) = h.α * x^2
conjugate_1d(h::L2norm, v::Float64) = v^2 / (4.0 * h.α)
prox_1d(h::L2norm, x::Float64, η::Float64) = x / (1.0 + 2.0 * η * h.α)
