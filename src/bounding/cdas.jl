"""
    CDAS

Coordinate Descent with Active-Set strategy solver for the lower and upper 
bounding steps in the [`BnbSolver`](@ref).

# Arguments

- `reltol::Float64` : Relative objective improvment tolearance for lower bouning.
- `maxiter_cd::Int` : Maximum number of cooridinate descent updates.
- `maxiter_as::Int` : Maximum number of active-set updates.
"""
struct CDAS <: AbstractBoundingSolver
    reltol::Float64
    maxiter_cd::Int
    maxiter_as::Int
    function CDAS(; 
        reltol::Float64=1e-4, 
        maxiter_cd::Int=1_000, 
        maxiter_as::Int=100
    )
        @assert reltol >= 0.
        @assert maxiter_cd >= 0
        @assert maxiter_as >= 0
        return new(reltol, maxiter_cd, maxiter_as)
    end
end

function cd_loop!(
    F::AbstractDatafit, 
    G::AbstractPenalty, 
    A::Matrix, 
    y::Vector, 
    x::Vector,
    w::Vector, 
    u::Vector, 
    ν::Vector, 
    ρ::Vector, 
    η::Vector, 
    δ::Vector, 
    S::BitArray, 
    Sb::BitArray,
    )
    for i in findall(S)
        ai = A[:, i]
        xi = x[i]
        ci = xi + ν[i] * (ai' * u)          # O(m)
        if Sb[i] & (abs(ci) <= δ[i])
            x[i] = prox_l1_1d(ci, η[i])     # O(1)
        else
            x[i] = prox_1d(G, ci, ρ[i])     # O(1)
        end
        if x[i] != xi
            axpy!(x[i] - xi, ai, w)         # O(m)
            copy!(u, -gradient(F, y, w))    # O(m)
        end
    end
    return nothing
end

function compute_primal_value(
    F::AbstractDatafit, 
    G::AbstractPenalty, 
    y::Vector, 
    λ::Float64,
    x::Vector,
    w::Vector, 
    τ::Float64,
    μ::Float64,
    S::BitArray,
    Sb::BitArray,
)
    Fval = value(F, y, w)                   # O(m)
    Gval = 0.
    for i in findall(S)
        xi = x[i]
        if Sb[i] & (abs(xi) <= μ)
            Gval += τ * abs(xi)             # O(1)
        else
            Gval += value_1d(G, xi) + 1.    # O(1)
        end
    end
    return Fval + λ * Gval
end

function compute_dual_value(
    F::AbstractDatafit, 
    G::AbstractPenalty, 
    A::Matrix,
    y::Vector, 
    λ::Float64,
    u::Vector,
    v::Vector,
    p::Vector,
    S::BitArray,
    Sb::BitArray;
)
    nz = findall(S)
    v[nz] = A[:, nz]' * u                               # O(m|nz|)
    p[nz] = conjugate_vectorized(G, v[nz] / λ) .- 1.
    cFval = conjugate(F, y, -u)
    cGval = 0.
    for i in nz
        cGval += Sb[i] ? max(p[i], 0.) : p[i]           # O(1)
    end
    return -cFval - λ * cGval
end

function update_active_set!(
    G::AbstractPenalty,
    A::Matrix,
    λ::Float64,
    u::Vector,
    v::Vector,
    p::Vector,
    τ::Float64,
    S::BitArray,
    Sbi::BitArray,
)
    # Only indices in Sbi may be added to S since indices of S1/Sbb ones are 
    # forced in S and indices of S0/Sb0 are excluded from S.
    violations = Vector{Int}()
    for i in findall(@. !S & Sbi)
        v[i] = A[:, i]' * u
        p[i] = conjugate_1d(G, v[i] / λ) .- 1.
        if abs(v[i]) > τ * λ
            push!(violations, i)
            S[i] = true
        end
    end
    return violations
end

function bound!(
    bounding_solver::CDAS, 
    problem::Problem, 
    solver::BnbSolver, 
    node::BnbNode, 
    options::BnbOptions,   
    bounding_type::BoundingType,
    )
    
    # ----- Initialization ----- #

    # Problem data
    F = problem.F
    G = problem.G
    A = problem.A
    y = problem.y
    λ = problem.λ
    a = problem.a
    m = problem.m
    n = problem.n

    # Node data
    if bounding_type == LOWER
        S0 = node.S0
        S1 = node.S1
        Sb = node.Sb
        x = node.x
        w = node.w
        u = node.u
    elseif bounding_type == UPPER
        S0 = copy(node.S0 .| node.Sb)
        S1 = copy(node.S1)
        Sb = falses(n)
        x = zeros(n)
        x[S1] = copy(node.x[S1])
        w = A[:, S1] * x[S1]
        u = -gradient(F, y, w)
    else
        error("Unknown bounding type $bounding_type")
    end

    S = (x .!= 0.)

    # Parameter values
    maxtime = options.maxtime
    tolgap = options.tolgap
    tolprune = options.tolprune
    dualpruning = options.dualpruning
    l0screening = options.l0screening
    l1screening = options.l1screening
    reltol = bounding_solver.reltol
    maxiter_cd = bounding_solver.maxiter_cd
    maxiter_as = bounding_solver.maxiter_as
   
    # Constants and working values
    τ = G.τ
    μ = G.μ
    α = lipschitz_constant(F, y)
    κ = α .* a
    ν = 1. ./ κ
    ρ = λ ./ κ
    η = τ .* ρ
    δ = η .+ μ
    v = Vector{Float64}(undef, n) 
    p = Vector{Float64}(undef, n) 

    # Active set and support configuration
    Sb0 = falses(n)
    Sbi = copy(Sb)
    Sbb = falses(n)
    S1i = copy(S1) 
    S1b = falses(n)
    # TODO : maybe keep S1 unsplitted

    # Objective values
    ub = solver.ub
    pv = Inf
    dv = -Inf

    # ----- Main loop ----- #

    it_as = 0
    it_cd = 0
    while true

        it_as += 1
        
        # ----- Inner CD solver ----- #

        while true
            it_cd += 1
            pv_old = pv
            cd_loop!(F, G, A, y, x, w, u, ν, ρ, η, δ, S, Sb)
            pv = compute_primal_value(F, G, y, λ, x, w, τ, μ, S, Sb)
            if bounding_type == LOWER
                (abs(pv - pv_old) / (abs(pv) + 1e-10) < reltol) && break
            elseif bounding_type == UPPER         
                dv = compute_dual_value(F, G, A, y, λ, u, v, p, S, Sb)
                (abs(pv - dv) / (abs(pv) + 1e-10) < tolgap) && break
            end
            (it_cd >= maxiter_cd) && break
        end

        # --- Active-set update and termination check --- #

        (bounding_type == UPPER) && break
        V = update_active_set!(G, A, λ, u, v, p, τ, S, Sbi) 
        (length(V) > 0) || break
        (it_as >= maxiter_as) && break
        (elapsed_time(solver) >= maxtime) && break

        # --- Accelerations --- #

        if (bounding_type == LOWER) & (l1screening | l0screening | dualpruning)
            dv = compute_dual_value(F, G, A, y, λ, u, v, p, S, Sb)
            gap = abs(pv - dv)
            l1screening && l1screening!(F, A, y, λ, x, w, u, v, α, τ, gap, Sb0, Sbi, Sbb)
            dualpruning && (dv >= ub + tolprune) && break
            l0screening && l0screening!(solver, node, F, A, y, λ, x, w, u, p, ub, dv, tolprune, S0, S1, Sb, Sbi, Sbb, S1i)
            S .|= S1
            S .|= Sbb
            S .&= (.!S0)
            S .&= (.!Sb0)
        end
    end

    # ----- Post-processing ----- #

    if bounding_type == LOWER
        node.lb = compute_dual_value(F, G, A, y, λ, u, v, p, S, Sb)
    else
        node.ub = pv
        node.x_ub = copy(x)
    end

    return nothing
end
