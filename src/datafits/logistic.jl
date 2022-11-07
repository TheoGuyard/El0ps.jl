struct Logistic <: AbstractDatafit end

Base.show(io::IO, F::Logistic) = print(io, "Logistic")
lipschitz_constant(F::Logistic, y::Vector) = 0.25 / length(y)

function value(F::Logistic, y::Vector, w::Vector)
    return sum(log.(1.0 .+ exp.(-y .* w))) / length(y)
end

function gradient(F::Logistic, y::Vector, w::Vector)
    return (@. -y / (1. + exp(y * w))) / length(y)
end

function conjugate(F::Logistic, y::Vector, u::Vector)
    v = -(u .* y) * length(y)
    r = 1. .- v
    return (v' * log.(v) + r' * log.(r)) / length(y)
end

function dual_scale!(F::Logistic, y::Vector, u::Vector)
    m = length(y)
    for j in eachindex(y, u)
        u[j] = clamp(u[j], min(y[j] / m, 0.), max(y[j] / m, 0.))
    end
    return nothing
end
