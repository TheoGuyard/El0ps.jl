"""
    CDAS

Coordinate Descent with Active-Set strategy solver for the lower and upper
bounding steps in the [`BnbSolver`](@ref).

# Arguments

- `bounding_type::BoundingType` : Bounding type.
- `reltol::Float64=1e-4` : Relative tolearance for lower bouning.
- `maxiter_cd::Int=1_000` : Maximum number of cooridinate descent updates.
- `maxiter_as::Int=100` : Maximum number of active-set updates.
"""
struct CDAS <: AbstractBoundingSolver
    bounding_type::BoundingType
    reltol::Float64
    maxiter_cd::Int
    maxiter_as::Int
    function CDAS(
        bounding_type::BoundingType;
        reltol::Float64 = 1e-4,
        maxiter_cd::Int = 1_000,
        maxiter_as::Int = 100,
    )
        @assert reltol >= 0.0
        @assert maxiter_cd >= 0
        @assert maxiter_as >= 0
        return new(bounding_type, reltol, maxiter_cd, maxiter_as)
    end
end

bounding_type(bounding_solver::CDAS) = bounding_solver.bounding_type

function cd_loop!(
    f::AbstractDatafit,
    h::AbstractPenalty,
    A::Matrix,
    x::Vector,
    w::Vector,
    u::Vector,
    ρ::Vector,
    η::Vector,
    δ::Vector,
    S::BitArray,
    Sb::BitArray,
)
    for i in findall(S)
        ai = A[:, i]
        xi = x[i]
        ci = xi + ρ[i] * (ai' * u)          # O(m)
        if Sb[i] & (abs(ci) <= δ[i])
            x[i] = prox_l1_1d(ci, η[i])     # O(1)
        else
            x[i] = prox_1d(h, ci, ρ[i])     # O(1)
        end
        if x[i] != xi
            axpy!(x[i] - xi, ai, w)         # O(m)
            copy!(u, -gradient(f, w))       # O(m)
        end
    end
    return nothing
end

function compute_primal_value(
    f::AbstractDatafit,
    h::AbstractPenalty,
    λ::Float64,
    x::Vector,
    w::Vector,
    τ::Float64,
    μ::Float64,
    S::BitArray,
    Sb::BitArray,
)
    fval = value(f, w)                      # O(m)
    hval = 0.0
    for i in findall(S)
        xi = x[i]
        if Sb[i] & (abs(xi) <= μ)
            hval += τ * abs(xi)             # O(1)
        else
            hval += value_1d(h, xi) + λ     # O(1)
        end
    end
    return fval + hval
end

function compute_dual_value(
    f::AbstractDatafit,
    h::AbstractPenalty,
    A::Matrix,
    λ::Float64,
    u::Vector,
    v::Vector,
    p::Vector,
    S::BitArray,
    Sb::BitArray;
)
    nz = findall(S)
    v[nz] = A[:, nz]' * u                               # O(m|nz|)
    p[nz] = [conjugate_1d(h, v[i]) .- λ for i in nz]
    cfval = conjugate(f, -u)
    chval = 0.0
    for i in nz
        chval += Sb[i] ? max(p[i], 0.0) : p[i]           # O(1)
    end
    return -cfval - chval
end

function update_active_set!(
    h::AbstractPenalty,
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
        p[i] = conjugate_1d(h, v[i]) .- λ
        if abs(v[i]) > τ
            push!(violations, i)
            S[i] = true
        end
    end
    return violations
end

function bound!(bounding_solver::CDAS, problem::Problem, solver, node, options)

    # ----- Initialization ----- #

    # Avoid recomputing the same upper bound two times
    if (bounding_solver.bounding_type == UPPER_BOUNDING) && (node.type == ZERO)
        node.ub = node.parent.ub
        node.x_ub = copy(node.parent.x_ub)
        node.status = SOLVED
        return nothing
    end

    # Problem data
    f = problem.f
    h = problem.h
    A = problem.A
    λ = problem.λ
    τ = problem.τ
    μ = problem.μ
    a = problem.a
    m = problem.m
    n = problem.n

    # Node data
    if bounding_solver.bounding_type == LOWER_BOUNDING
        S0 = node.S0
        S1 = node.S1
        Sb = node.Sb
        x = node.x
        w = node.w
        u = node.u
    elseif bounding_solver.bounding_type == UPPER_BOUNDING
        S0 = copy(node.S0 .| node.Sb)
        S1 = copy(node.S1)
        Sb = falses(n)
        x = zeros(n)
        x[S1] = copy(node.x[S1])
        w = A[:, S1] * x[S1]
        u = -gradient(f, w)
    else
        error("Unknown bounding type $bounding_type")
    end

    # Parameter values
    maxtime = options.maxtime
    tolgap = options.tolgap
    tolprune = options.tolprune
    dualpruning = options.dualpruning
    l0screening = options.l0screening
    l1screening = options.l1screening
    reltol = bounding_solver.reltol
    cdltol = 0.5 * reltol
    maxiter_cd = bounding_solver.maxiter_cd
    maxiter_as = bounding_solver.maxiter_as

    # Constants and working values
    α = lipschitz_constant(f)
    κ = α .* a
    ρ = 1.0 ./ κ
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
    S = @. (x != 0.0) | S1
    # TODO : maybe keep S1 unsplitted

    # Objective values
    ub = solver.ub
    pv = Inf
    dv = NaN

    # ----- Main loop ----- #

    it_as = 0
    it_cd = 0
    while true

        it_as += 1
        v .= NaN
        p .= NaN
        dv = NaN

        # ----- Inner CD solver ----- #

        while true
            it_cd += 1
            pv_old = pv
            cd_loop!(f, h, A, x, w, u, ρ, η, δ, S, Sb)
            pv = compute_primal_value(f, h, λ, x, w, τ, μ, S, Sb)
            if bounding_solver.bounding_type == LOWER_BOUNDING
                (abs(pv - pv_old) / (abs(pv) + 1e-10) < cdltol) && break
            elseif bounding_solver.bounding_type == UPPER_BOUNDING
                dv = compute_dual_value(f, h, A, λ, u, v, p, S, Sb)
                (abs(pv - dv) / (abs(pv) + 1e-10) < tolgap) && break
            end
            (it_cd >= maxiter_cd) && break
        end

        # --- Active-set update and termination check --- #

        (bounding_solver.bounding_type == UPPER_BOUNDING) && break
        V = update_active_set!(h, A, λ, u, v, p, τ, S, Sbi)
        (it_as >= maxiter_as) && break
        (elapsed_time(solver) >= maxtime) && break
        if length(V) == 0
            dv = compute_dual_value(f, h, A, λ, u, v, p, S, Sb)
            (abs(pv - dv) / (abs(pv) + 1e-10) < reltol) && break
            (cdltol <= 1e-8) && break
            cdltol *= 1e-2
        end

        # --- Accelerations --- #

        if (bounding_solver.bounding_type == LOWER_BOUNDING) &
           (l1screening | l0screening | dualpruning)
            dv = isnan(dv) ? compute_dual_value(f, h, A, λ, u, v, p, S, Sb) : dv
            gap = abs(pv - dv)
            dualpruning && (dv >= ub + tolprune) && break
            l1screening && l1screening!(f, A, λ, x, w, u, v, α, τ, gap, Sb0, Sbi, Sbb)
            l0screening && l0screening!(
                f,
                A,
                λ,
                x,
                w,
                u,
                p,
                ub,
                dv,
                tolprune,
                S0,
                S1,
                Sb,
                Sbi,
                Sbb,
                S1i,
            )
            S .|= S1
            S .|= Sbb
            S .&= (.!S0)
            S .&= (.!Sb0)
        end
    end

    # ----- Post-processing ----- #

    if bounding_solver.bounding_type == LOWER_BOUNDING
        node.lb = isnan(dv) ? compute_dual_value(f, h, A, λ, u, v, p, S, Sb) : dv
        node.lb_it = it_cd
        node.lb_l1screening_Sb0 = sum(Sb0)
        node.lb_l1screening_Sbb = sum(Sbb)
        node.lb_l0screening_S0 =
            sum(S0) - (isa(node.parent, Nothing) ? 1 : (sum(node.parent.S0) + 1))
        node.lb_l0screening_S1 =
            sum(S1) - (isa(node.parent, Nothing) ? 1 : (sum(node.parent.S1) + 1))
    elseif bounding_solver.bounding_type == UPPER_BOUNDING
        node.ub = pv
        node.x_ub = copy(x)
        node.status = SOLVED
    end

    return nothing
end
