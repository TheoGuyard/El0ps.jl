"""
    BigmL2norm <: AbstractPenalty

Big-M constraint plus L2-norm function
`h(x) = h.α * norm(x, 2)^2 + (norm(x, Inf) < h.M ? 0. : Inf)`, where `h.α > 0`
and `h.M > 0`.

# Arguments

- `M::Float64` : Big-M value.
- `α::Float64` : L2 regularization strength.
"""
struct BigmL2norm <: AbstractPenalty
    M::Float64
    α::Float64
    function BigmL2norm(M::Float64, α::Float64)
        (M > 0.0) || error("Parameter M must be positive")
        (α > 0.0) || error("Parameter α must be positive")
        return new(M, α)
    end
end

Base.show(io::IO, h::BigmL2norm) = print(io, "Bigm + L2-norm")
compute_τ(h::BigmL2norm) =
    (sqrt(1.0 / h.α) < h.M) ? sqrt(4.0 * h.α) : (1.0 / h.M) + h.α * h.M
compute_μ(h::BigmL2norm) = (sqrt(1.0 / h.α) < h.M) ? sqrt(1.0 / h.α) : h.M
value_1d(h::BigmL2norm, x::Float64) = abs(x) <= h.M ? h.α * x^2 : Inf
function conjugate_1d(h::BigmL2norm, v::Float64)
    r = clamp(v / (2.0 * h.α), -h.M, h.M)
    return v * r - h.α * r^2
end
prox_1d(h::BigmL2norm, x::Float64, η::Float64) = clamp(x / (1.0 + 2.0 * η * h.α), -h.M, h.M)
dual_scaling_factor(h::BigmL2norm, v::Vector) = 1.0
