"""
    L2norm <: AbstractPerturbation

L2-norm function `h(x) = h.α * norm(x, 2)^2`, where `h.α > 0`.

# Arguments

- `α::Float64` : L2 regularization strength.
"""
struct L2norm <: AbstractPerturbation
    α::Float64
    τ::Float64
    μ::Float64
end

"""
    L2norm(α::Float64)

[`L2norm`](@ref) constructor.
"""
function L2norm(α::Float64)
    (α > 0.0) || error("Parameter α must be positive")
    τ = 2.0 * sqrt(α)
    μ = 1.0 / sqrt(α)
    return L2norm(α, τ, μ)
end

Base.show(io::IO, h::L2norm) = print(io, "L2-norm")
value_1d(h::L2norm, x::Float64) = h.α * x^2
conjugate_1d(h::L2norm, v::Float64) = v^2 / (4.0 * h.α)
prox_1d(h::L2norm, x::Float64, η::Float64) = x / (1.0 + 2.0 * η * h.α)
