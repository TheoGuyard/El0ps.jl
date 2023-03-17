"""
    NegLogSymtri <: AbstractPenalty

Function `h(x)= -α * ln(1-|x|/σ) + I(|x| ≤ σ)`, which is `α` times the
Negative-Loglikelihood of the symmetric triangular distributions.

# Arguments

- `α::Float64` : Multiplicative factor.
- `σ::Float64` : Distribution spread parameter.
"""
struct NegLogSymtri <: AbstractPenalty
    α::Float64
    σ::Float64
    function NegLogSymtri(α::Float64, σ::Float64)
        (α > 0.0) || error("Parameter α must be positive")
        (σ > 0.0) || error("Parameter σ must be positive")
        return new(α, σ)
    end
end

Base.show(io::IO, h::NegLogSymtri) = print(io, "Neg-Log of sym-tri. distrib.")
compute_τ(h::NegLogSymtri) = approximate_τ(h)
compute_μ(h::NegLogSymtri) = (σ / α) - 1.0 / compute_τ(h)
value_1d(h::NegLogSymtri, x::Float64) = abs(x) <= h.σ ? -h.α * log(1. - abs(x) / h.σ) : Inf
function conjugate_1d(h::NegLogSymtri, v::Float64)
    u = max((σ / α) * abs(v) - 1.0, 0.0)
    return h.α * (u - log(u + 1.0))
end
function prox_1d(h::NegLogSymtri, x::Float64, η::Float64)
    a = -x - h.σ * sign(x)
    b = h.σ * abs(x) - η * h.α
    z = 0.5 * (sign(x) * sqrt((h.σ - abs(x))^2 + 4.0 * η * h.α) - a)
    return clamp(z, -h.σ, h.σ)
end
