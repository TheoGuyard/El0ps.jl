function sparse_vector(n::Int, k::Int)
    x = zeros(Float64, n)
    s = Int.(floor.(collect(1:n/k:n)))
    x[s] = sign.(rand(Normal(0, 1.0), k))
    return x
end

function correlated_dictionary(
    m::Int, 
    n::Int, 
    ρ::Float64; 
    normalize::Bool=false
)
    μ = zeros(n)
    N = [collect(1:n);]
    Σ = ρ .^ (abs.(repeat(N, inner=(1, n)) - repeat(N', inner=(n, 1))))
    L = MvNormal(μ, Σ)
    A = zeros(m, n)
    for j in 1:m
        A[j, :] = rand(L)
    end
    if normalize
        for ai in eachcol(A)
            normalize!(ai)
        end
    end
    return A
end

"""
    synthetic_data_regression(
        k::Int,
        m::Int,
        n::Int,
        ρ::Float64,
        σ::Float64;
        normalize::Bool=false
    )

Generate a synthetic regression dataset `y = Ax + ϵ` where `x` is a k-sparse 
vector.

# Arguments

- `k::Int` : Number of non-zeros in `x`.
- `m::Int` : Number of rows in `A`.
- `n::Int` : Number of columns in `A`.
- `ρ::Float64` : Correlation parameter `ρ ∈ [0,1]` between the columns in `A`.
- `σ::Float64` : SNR ratio of `ϵ` with respect to `Ax`.
- `normalize::Bool=false` : Whether to normalize the columns in `A`.
"""
function synthetic_data_regression(
    k::Int,
    m::Int,
    n::Int,
    ρ::Float64,
    σ::Float64;
    normalize::Bool=false
)

    # Consistence checks
    (0 <= k <= n) || throw(ArgumentError("Parameter k must be between 0 and n"))
    (0.0 <= ρ <= 1.0) || throw(ArgumentError("Parameter ρ must be between 0 and 1"))
    (0.0 <= σ) || throw(ArgumentError("Parameter σ must be positive"))

    x = sparse_vector(n, k)
    A = correlated_dictionary(m, n, ρ, normalize=normalize)
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
        normalize::Bool=false
    )

Generate a synthetic classifiaction dataset where `yj = 1` with proba 
`1 / (1 + exp(-σ * [Ax]_j))` and `yj = -1` otheriwse.

# Arguments

- `k::Int` : Number of non-zeros in `x`.
- `m::Int` : Number of rows in `A`.
- `n::Int` : Number of columns in `A`.
- `ρ::Float64` : Correlation parameter `ρ ∈ [0,1]` between the columns in `A`.
- `σ::Float64` : Probability parameter.
- `normalize::Bool=false` : Whether to normalize the columns in `A`.
"""
function synthetic_data_classification(
    k::Int,
    m::Int,
    n::Int,
    ρ::Float64,
    σ::Float64;
    normalize::Bool=false
)

    # Consistence checks
    (0 <= k <= n) || throw(ArgumentError("Parameter k must be between 0 and n"))
    (0.0 <= ρ <= 1.0) || throw(ArgumentError("Parameter ρ must be between 0 and 1"))
    (0.0 <= σ) || throw(ArgumentError("Parameter σ must be positive"))

    x = sparse_vector(n, k)
    A = correlated_dictionary(m, n, ρ, normalize=normalize)
    w = A * x
    p = @. 1.0 / (1.0 + exp(-σ * w))
    y = @. 2.0 * Float64(rand(Bernoulli(p))) - 1.0
    
    return x, A, y
end