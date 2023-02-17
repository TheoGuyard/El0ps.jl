"""
    L1L2norm

L1L2-norm function h(x) = α * ||x||_1 + β * ||x||_2^2.

# Arguments

- `α::Float64` : L1 regularization strength.
- `β::Float64` : L2 regularization strength.
"""
struct L1L2norm <: AbstractPerturbation
    α::Float64
    β::Float64
    τ::Float64
    μ::Float64
    function L1L2norm(α::Float64, β::Float64)
        (α > 0.) || error("Parameter α must be positive")
        (β > 0.) || error("Parameter β must be positive")
        τ =  α + sqrt(4. * β)
        μ = sqrt(1. / β)
        return new(α, β, τ, μ)
    end
end

Base.show(io::IO, h::L1L2norm) = print(io, "L1L2-norm")
value_1d(h::L1L2norm, x::Float64) = h.α * abs(x) + h.β * x^2
conjugate_1d(h::L1L2norm, v::Float64) = max(abs(v) - h.α, 0.)^2 / (4. * h.β)
prox_1d(h::L1L2norm, x::Float64, η::Float64) = (sign(x) / (1. + 2. * η * h.β)) * max(abs(x) - η * h.α, 0.)