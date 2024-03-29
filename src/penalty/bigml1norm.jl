"""
    BigmL1norm <: AbstractPenalty

Big-M constraint plus L1-norm function
`h(x) = h.α * norm(x, 1) + (norm(x, Inf) < h.M ? 0. : Inf)`, where `h.α > 0` and
`h.M > 0`.

# Arguments

- `M::Float64` : Big-M value.
- `α::Float64` : L1 regularization strength.
"""
struct BigmL1norm <: AbstractPenalty
    M::Float64
    α::Float64
    function BigmL1norm(M::Float64, α::Float64)
        (M > 0.0) || error("Parameter M must be positive")
        (α > 0.0) || error("Parameter α must be positive")
        return new(M, α)
    end
end

Base.show(io::IO, h::BigmL1norm) = print(io, "Bigm + L1-norm")
compute_τ(h::BigmL1norm) = (1.0 / h.M) + h.α
compute_μ(h::BigmL1norm) = h.M
value_1d(h::BigmL1norm, x::Float64) = abs(x) <= h.M ? h.α * abs(x) : Inf
conjugate_1d(h::BigmL1norm, v::Float64) = h.M * max(abs(v) - h.α, 0.0)
prox_1d(h::BigmL1norm, x::Float64, η::Float64) = sign(x) * clamp(abs(x) - η * h.α, 0.0, h.M)
dual_scaling_factor(h::BigmL1norm, v::Vector) = 1.0
