"""
    CD

Coordinate Descent solver for the lower and upper bounding steps in the 
[`BnbSolver`](@ref).

# Arguments

- `tolgap::Float64` : Absolute tolearance on the duality gap.
- `maxiter::Int` : Maximum number of itrations.
"""
struct CD <: AbstractBoundingSolver
    tolgap::Float64
    maxiter::Int
    function CD(; tolgap::Float64=1e-8, maxiter::Int=10_000)
        @assert tolgap >= 0.
        @assert maxiter >= 0
        return new(tolgap, maxiter)
    end
end

function bound!(
    bounding_solver::CD, 
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

    # Recover node informations
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
    

    # Support configuration
    Sb0 = falses(n)
    Sbi = copy(Sb)
    Sbb = falses(n)
    S1i = copy(S1)
    S1b = falses(n)

    # BnB, options and solver values
    ub = solver.ub
    tolprune = options.tolprune
    maxtime = options.maxtime
    dualpruning = options.dualpruning
    l0screening = options.l0screening
    l1screening = options.l1screening
    tolgap = bounding_solver.tolgap
    maxiter = bounding_solver.maxiter

    # Useful constants
    τ = G.τ
    μ = G.μ
    α = lipschitz_constant(F, y)
    κ = α .* a
    ν = 1. ./ κ
    ρ = λ ./ κ
    η = τ .* ρ
    δ = η .+ μ
    

    # Objectives
    pv  = Inf
    dv  = -Inf
    gap = Inf


    # ----- Main loop ----- #

    it = 0
    while true

        it += 1

        # ----- Descent loop ----- #

        idx = (Sbi .| Sbb .| S1i .| S1b)
        for i in shuffle(findall(idx))
            ai = A[:, i]
            xi = x[i]
            ci = xi + ν[i] * (ai' * u)
            if Sbi[i] | Sbb[i]
                if abs(ci) < δ[i]
                    x[i] = prox_l1_1d(ci, η[i])
                else
                    x[i] = prox_1d(G, ci, ρ[i])
                end
            else
                x[i] = prox_1d(G, ci, ρ[i])
            end
            if x[i] != xi
                axpy!(x[i] - xi, ai, w)
                copy!(u, -gradient(F, y, w))
            end
        end


        # ----- Gap computation ----- #

        dual_scale!(F, y, u)
        v = fill(NaN, n)
        p = fill(NaN, n)
        v[idx] = dual_scale!(G, A[:, idx], u, λ)
        p[idx] = conjugate_vectorized(G, v[idx] ./ λ) .- 1.

        xSb       = x[Sbi .| Sbb]
        xSb_below = xSb[abs.(xSb) .< μ]
        xSb_above = xSb[abs.(xSb) .>= μ]
        μSbi      = any(Sbi) ? sum(max.(p[Sbi], 0.)) : 0.
        μSbb      = any(Sbb) ? sum(max.(p[Sbb], 0.)) : 0.
        μS1       = any(S1) ? sum(p[S1]) : 0.

        pv = value(F, y, w) + λ * (
            τ * norm(xSb_below, 1) +
            value(G, xSb_above) + length(xSb_above) +
            value(G, x[S1]) + sum(S1)
        )
        dv = -conjugate(F, y, -u) - λ * (μSbi + μSbb + μS1)
        gap = abs(pv - dv)

        # ----- Stopping criterion ----- #

        (gap / (pv + 1e-10) < tolgap) && break
        (it > maxiter) && break
        (elapsed_time(solver) >= maxtime) && break
        
        # --- Accelerations --- #
        
        if bounding_type == LOWER
            dualpruning && (dv >= ub + tolprune) && break
            l0screening && l0screening!(solver, node, F, A, y, λ, x, w, u, p, ub, dv, tolprune, S0, S1, Sb, Sbi, Sbb, S1i)
            l1screening && l1screening!(F, A, y, λ, x, w, u, v, α, τ, gap, Sb0, Sbi, Sbb)
        end
    end

    # ----- Post-processing ----- #

    if bounding_type == LOWER
        node.lb = dv
        node.lb_it = it
        node.lb_l1screening_Sb0 = sum(Sb0)
        node.lb_l1screening_Sbb = sum(Sbb)
        node.lb_l0screening_S0 = sum(S0) - (isa(node.parent, Nothing) ? 1 : (sum(node.parent.S0) + 1))
        node.lb_l0screening_S1 = sum(S1) - (isa(node.parent, Nothing) ? 1 : (sum(node.parent.S1) + 1))
    elseif bounding_type == UPPER
        node.ub = pv
        node.x_ub = copy(x)
    else
        error("Unknown bounding type $bounding_type")
    end

    return nothing
end
