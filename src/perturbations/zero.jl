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

Base.show(io::IO, h::Zero) = print(io, "Zero-function perturbation")
value_1d(h::Zero, x::Float64) = 0.
conjugate_1d(h::Zero, v::Float64) = (v == 0.) ? 0. : +Inf
prox_1d(h::Zero, x::Float64, η::Float64) = x
