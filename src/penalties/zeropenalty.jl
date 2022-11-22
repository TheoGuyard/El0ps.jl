"""
    ZeroPenalty <: AbstractPenalty 

Zero-function

``G(\\mathbf{x}) = 0``

This penalty is to be used with a [`DirectSolver`](@ref) only in order to force 
the use of a SOS formulation.
"""
struct ZeroPenalty <: AbstractPenalty 
    τ::Float64
    μ::Float64
    ZeroPenalty() = new(Inf, Inf)
end

Base.show(io::IO, G::ZeroPenalty) = print(io, "Zero penalty")

value_1d(G::ZeroPenalty, x::Float64) = 0.
conjugate_1d(G::ZeroPenalty, v::Float64) = (v == 0.) ? 0. : +Inf
prox_1d(G::ZeroPenalty, x::Float64, η::Float64) = x
dual_scale!(G::ZeroPenalty, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(G::ZeroPenalty, model::JuMP.Model)
    for i in eachindex(model[:x], model[:z])
        @constraint(model, [model[:x][i], model[:z][i]] in SOS1())
    end
    @constraint(model, model[:Ωcost] >= sum(model[:z]))
end
