struct L1L2norm <: AbstractPenalty
    α::Float64
    β::Float64
    τ::Float64
    μ::Float64
    function L1L2norm(α::Float64, β::Float64)
        (α > 0.) || error("Parameter α must be positive")
        (β > 0.) || error("Parameter β must be positive")
        τ =  α + sqrt(2. * β)
        μ = sqrt(2. / β)
        return new(α, β, τ, μ)
    end
end

Base.show(io::IO, G::L1L2norm) = print(io, "L1L2-norm")
value_1d(G::L1L2norm, x::Float64) = G.α * abs(x) + (G.β / 2.) * x^2
conjugate_1d(G::L1L2norm, v::Float64) = max(abs(v) - G.α, 0.)^2 / (2. * G.β)
prox_1d(G::L1L2norm, x::Float64, η::Float64) = (sign(x) / (1. + η * G.β)) * max(abs(x) - η * G.α, 0.)
dual_scale!(G::L1L2norm, A::Matrix, u::Vector, λ::Float64) = A' * u
