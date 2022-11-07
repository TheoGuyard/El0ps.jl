struct BigmL1norm <: AbstractPenalty
    M::Float64
    α::Float64
    τ::Float64
    μ::Float64
    function BigmL1norm(M::Float64, α::Float64)
        (M > 0.) || error("Parameter M must be positive")
        (α > 0.) || error("Parameter α must be positive")
        return new(M, α, (1. / M) + α, M)
    end
end

Base.show(io::IO, G::BigmL1norm) = print(io, "Bigm + L1-norm")
value_1d(G::BigmL1norm, x::Float64) = abs(x) <= G.M ? G.α * abs(x) : Inf
conjugate_1d(G::BigmL1norm, v::Float64) = G.M * (max(abs(v) - G.α, 0.))
prox_1d(G::BigmL1norm, x::Float64, η::Float64) = sign(x) * clamp(abs(x) - η * G.α, 0., G.M)
dual_scale!(G::BigmL1norm, A::Matrix, u::Vector, λ::Float64) = A' * u
