struct L1norm <: AbstractPenalty
    α::Float64
    τ::Float64
    μ::Float64
    function L1norm(α::Float64)
        (α > 0.) || error("Parameter α must be positive")
        return new(α, α, Inf)
    end
end

Base.show(io::IO, G::L1norm) = print(io, "L1-norm")
value_1d(G::L1norm, x::Float64) = G.α * abs(x)
conjugate_1d(G::L1norm, v::Float64) = (abs.(v) <= G.α) ? 0. : Inf
prox_1d(G::L1norm, x::Float64, η::Float64) = sign(x) * max(abs(x) - η * G.α, 0.0)
function dual_scale!(G::L1norm, A::Matrix, u::Vector, λ::Float64)
    v = A' * u
    s = λ * G.α / (norm(v, Inf) + 1e-16)
    u *= s
    v *= s
    return v
end

function bind_model!(G::L1norm, model::JuMP.Model)
    n = length(model[:x])
    @variable(model, xabs[1:n])
    @constraint(model, xabs .* model[:z] .>= model[:x])
    @constraint(model, xabs .* model[:z] .>= -model[:x])
    @constraint(model, model[:Ωcost] >= sum(model[:z]) + G.α * sum(xabs))
    return nothing
end

