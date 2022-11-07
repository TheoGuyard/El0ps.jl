struct L2norm <: AbstractPenalty
    α::Float64
    τ::Float64
    μ::Float64
    function L2norm(α::Float64)
        (α > 0.) || error("Parameter α must be positive")
        return new(α, 2. * sqrt(α), 1. / sqrt(α))
    end
end

Base.show(io::IO, G::L2norm) = print(io, "L2-norm")
value_1d(G::L2norm, x::Float64) = G.α * x^2
conjugate_1d(G::L2norm, v::Float64) = v^2 / (4. * G.α)
prox_1d(G::L2norm, x::Float64, η::Float64) = x / (1. + 2. * η * G.α)
dual_scale!(G::L2norm, A::Matrix, u::Vector, λ::Float64) = A' * u
