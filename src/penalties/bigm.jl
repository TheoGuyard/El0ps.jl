"""
    Bigm <: AbstractPerturbation 

Convex indicator of a Big-M constraint

``h(\\mathbf{x}) = \\mathrm{Ind}(\\|\\mathbf{x}\\|_{\\infty} \\leq M)``

where `M > 0`.

# Arguments

- `M::Float64` : Big-M value.
"""
struct Bigm <: AbstractPerturbation 
    M::Float64
    τ::Float64
    μ::Float64
    Bigm(M::Float64) = new(M, 1. / M, M)
end

Base.show(io::IO, h::Bigm) = print(io, "Big-M constraint")
value_1d(h::Bigm, x::Float64) = abs(x) <= h.M ? 0. : Inf
conjugate_1d(h::Bigm, v::Float64) = h.M * abs(v)
prox_1d(h::Bigm, x::Float64, η::Float64) = clamp(x, -h.M, h.M)
dual_scale!(h::Bigm, A::Matrix, u::Vector, λ::Float64) = A' * u

function bind_model!(h::Bigm, model::JuMP.Model)
    @constraint(model, model[:x] .<= h.M .* model[:z])
    @constraint(model, model[:x] .>= -h.M .* model[:z])
    @constraint(model, model[:Ωcost] >= sum(model[:z]))
end
