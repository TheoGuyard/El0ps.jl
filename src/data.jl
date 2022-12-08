function sparse_vector(n::Int, k::Int)
    x = zeros(Float64, n)
    s = Int.(floor.(collect(1:n/k:n)))
    x[s] = sign.(rand(Normal(0, 1.0), k))
    return x
end

function correlated_dictionary(m::Int, n::Int, ρ::Float64)
    μ = zeros(n)
    N = [collect(1:n);]
    Σ = ρ .^ (abs.(repeat(N, inner=(1, n)) - repeat(N', inner=(n, 1))))
    L = MvNormal(μ, Σ)
    A = zeros(m, n)
    for j in 1:m
        A[j, :] = rand(L)
    end
    return A
end

function normalize!(A::Matrix)
    for ai in eachcol(A)
        ai .-= mean(ai)
        ai ./= norm(ai)
    end
    return nothing
end

"""
    synthetic_data_regression(
        k::Int,
        m::Int,
        n::Int,
        ρ::Float64,
        σ::Float64;
    )

Generate a synthetic regression instance of a [`Problem`](@ref) with 
`y = Ax + ϵ` where `x` is a k-sparse vector.

# Arguments

- `k::Int` : Number of non-zeros in `x`.
- `m::Int` : Number of rows in `A`.
- `n::Int` : Number of columns in `A`.
- `ρ::Float64` : Correlation between the columns in `A`.
- `σ::Float64` : SNR ratio of `ϵ` with respect to `Ax`.
"""
function synthetic_data_regression(
    k::Int,
    m::Int,
    n::Int,
    ρ::Float64,
    σ::Float64;
)

    @assert 0 <= k <= n
    @assert 0.0 <= ρ <= 1.0
    @assert 0.0 <= σ

    x = sparse_vector(n, k)
    A = correlated_dictionary(m, n, ρ)
    normalize!(A)
    w = A * x
    ϵ = rand(Normal(), m)
    ϵ *= ((σ != Inf) ? norm(w, 2) / (sqrt(σ) * norm(ϵ, 2)) : 0.0)
    y = w + ϵ

    return x, A, y
end

"""
    synthetic_data_classification(
        k::Int,
        m::Int,
        n::Int,
        ρ::Float64,
        σ::Float64;
    )

Generate a synthetic classification instance of a [`Problem`](@ref) with 
`y ≈ Bernoulli(p)` where `y[j] = 1` with probability 
`p[j] = 1 / (1 + exp(-σ[Ax]_j))` and `y[j] = -1` otherwise.

# Arguments

- `k::Int` : Number of non-zeros in `x`.
- `m::Int` : Number of rows in `A`.
- `n::Int` : Number of columns in `A`.
- `ρ::Float64` : Correlation between the columns in `A`.
- `σ::Float64` : Probability parameter.
"""
function synthetic_data_classification(
    k::Int,
    m::Int,
    n::Int,
    ρ::Float64,
    σ::Float64;
)

    @assert 0 <= k <= n
    @assert 0.0 <= ρ <= 1.0
    @assert 0.0 <= σ

    x = sparse_vector(n, k)
    A = correlated_dictionary(m, n, ρ)
    w = A * x
    p = @. 1.0 / (1.0 + exp(-σ * w))
    y = @. 2.0 * Float64(rand(Bernoulli(p))) - 1.0
    
    return x, A, y
end
