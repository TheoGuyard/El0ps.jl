"""
    L1norm <: AbstractPerturbation

L1-norm function `h(x) = h.α * norm(x, 1)`, where `h.α > 0`.

# Arguments

- `α::Float64` : L1 regularization strength.
"""
struct L1norm <: AbstractPerturbation
    α::Float64
    τ::Float64
    μ::Float64
end

"""
    L1norm(α::Float64)

[`L1norm`](@ref) constructor.
"""
function L1norm(α::Float64)
    (α > 0.0) || error("Parameter α must be positive")
    return L1norm(α, α, Inf)
end

Base.show(io::IO, h::L1norm) = print(io, "L1-norm")
value_1d(h::L1norm, x::Float64) = h.α * abs(x)
conjugate_1d(h::L1norm, v::Float64) = (abs.(v) <= h.α) ? 0.0 : Inf
prox_1d(h::L1norm, x::Float64, η::Float64) = sign(x) * max(abs(x) - η * h.α, 0.0)
