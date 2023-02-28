"""
    Bigm <: AbstractPerturbation

Convex indicator of a Big-M constraint `abs.(x) .<= h.M`, where `h.M > 0`.

# Arguments

- `M::Float64` : Big-M value.
"""
struct Bigm <: AbstractPerturbation
    M::Float64
    function Bigm(M::Float64)
        (M > 0.0) || error("Parameter M must be positive")
        return new(M)
    end
end

Base.show(io::IO, h::Bigm) = print(io, "Big-M constraint")
compute_τ(h::Bigm, λ::Float64) = λ / h.M
compute_μ(h::Bigm, λ::Float64) = h.M
value_1d(h::Bigm, x::Float64) = abs(x) <= h.M ? 0.0 : Inf
conjugate_1d(h::Bigm, v::Float64) = h.M * abs(v)
prox_1d(h::Bigm, x::Float64, η::Float64) = clamp(x, -h.M, h.M)
