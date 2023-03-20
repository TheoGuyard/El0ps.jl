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
compute_μ(h::NegLogSymtri) = h.σ - 1.0 / compute_τ(h)
value_1d(h::NegLogSymtri, x::Float64) = abs(x) <= h.σ ? -h.α * log(1. - abs(x) / h.σ) : Inf
function conjugate_1d(h::NegLogSymtri, v::Float64)
    u = max((h.σ / h.α) * abs(v) - 1.0, 0.0)
    return h.α * (u - log(u + 1.0))
end
function prox_1d(h::NegLogSymtri, x::Float64, η::Float64)
    z = h.σ - sqrt((h.σ - abs(x))^2 + 4.0 * η * h.α)
    return clamp(0.5 * (x + sign(x) * z), -h.σ, h.σ)
end
