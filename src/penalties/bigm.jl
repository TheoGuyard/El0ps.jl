"""
    Bigm

Big-M function G(x) = Ind(||x||_Inf <= M).

# Arguments

- `M::Float64` : Big-M value.
"""
struct Bigm <: AbstractPenalty 
    M::Float64
    τ::Float64
    μ::Float64
    Bigm(M::Float64) = new(M, 1. / M, M)
end

Base.show(io::IO, G::Bigm) = print(io, "Big-M constraint")
value_1d(G::Bigm, x::Float64) = abs(x) <= G.M ? 0. : Inf
conjugate_1d(G::Bigm, v::Float64) = G.M * abs(v)
prox_1d(G::Bigm, x::Float64, η::Float64) = clamp(x, -G.M, G.M)
dual_scale!(G::Bigm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(G::Bigm, model::JuMP.Model)
    @constraint(model, model[:x] .<= G.M .* model[:z])
    @constraint(model, model[:x] .>= -G.M .* model[:z])
    @constraint(model, model[:Ωcost] >= sum(model[:z]))
end