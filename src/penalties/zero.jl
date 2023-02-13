"""
    Zero <: AbstractPerturbation 

Zero-function

``h(\\mathbf{x}) = 0``

This perturbation is to be used with a [`DirectSolver`](@ref) only in order to force 
the use of a SOS formulation.
"""
struct Zero <: AbstractPerturbation 
    τ::Float64
    μ::Float64
    Zero() = new(Inf, Inf)
end

Base.show(io::IO, h::Zero) = print(io, "Zero perturbation")

value_1d(h::Zero, x::Float64) = 0.
conjugate_1d(h::Zero, v::Float64) = (v == 0.) ? 0. : +Inf
prox_1d(h::Zero, x::Float64, η::Float64) = x

function dual_scale!(h::Zero, A::Matrix, u::Vector, λ::Float64)
    copy!(u, zeros(size(u)))
    return A' * u
end

function bind_model!(h::Zero, model::JuMP.Model)
    for i in eachindex(model[:x], model[:z])
        @constraint(model, [model[:x][i], model[:z][i]] in SOS1())
    end
    @constraint(model, model[:Gcost] >= sum(model[:z]))
end
