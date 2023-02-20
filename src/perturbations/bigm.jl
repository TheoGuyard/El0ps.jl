"""
    Bigm <: AbstractPerturbation 

Convex indicator of a Big-M constraint `abs.(x) .<= h.M`, where `h.M > 0`.

# Arguments

- `M::Float64` : Big-M value.
"""
struct Bigm <: AbstractPerturbation 
    M::Float64
    τ::Float64
    μ::Float64
end

"""
    Bigm(M::Float64)

[`Bigm`](@ref) constructor.
"""
Bigm(M::Float64) = Bigm(M, 1. / M, M)

Base.show(io::IO, h::Bigm) = print(io, "Big-M constraint")
value_1d(h::Bigm, x::Float64) = abs(x) <= h.M ? 0. : Inf
conjugate_1d(h::Bigm, v::Float64) = h.M * abs(v)
prox_1d(h::Bigm, x::Float64, η::Float64) = clamp(x, -h.M, h.M)
