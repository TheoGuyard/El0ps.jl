"""
    ExplorationStrategy

Exploration strategy of a [`BnbSolver`](@ref).

- `BFS` : Breadth-First Search
- `DFS` : Depth-First Search
"""
@enum ExplorationStrategy begin
    BFS
    DFS
end

"""
    BranchingStrategy

Branching strategy of a [`BnbSolver`](@ref).

- LARGEST : Select the largest index in absolute value in the relaxed solution.
"""
@enum BranchingStrategy begin
    LARGEST
end

"""
    BnbOptions

Options of a [`BnbSolver`](@ref).

# Arguments 

- `lb_solver::AbstractBoundingSolver` : Solver for the lower-bounding step.
- `ub_solver::AbstractBoundingSolver` : Solver for the upper-bounding step.
- `exploration::ExplorationStrategy` : Tree exploration strategy.
- `branching::BranchingStrategy` : Tree branching strategy.
- `maxtime::Float64` : Maximum solution time in seconds.
- `maxnode::Int` : Maximum number of nodes.
- `tolgap::Float64` : Relative MIP gap tolerance.
- `tolint::Float64` : Integer tolerance, i.e., `x = 0` when `|x| < tolint`.
- `tolprune::Float64` : Prune a `node` in the `bnb` tree when `bnb.ub + tolprune < node.lb`.
- `dualpruning::Bool` : Toogle the dual-pruning acceleration.
- `l0screening::Bool` : Toogle the L0-screening acceleration.
- `l1screening::Bool` : Toogle the L1-screening acceleration.
- `verbosity::Bool` : Toogle verbosity.
- `showevery::Int` : Displays logs every `showevery` nodes explored.
- `keeptrace::Bool` : Whether to fill the [`BnbTrace`](@ref) or not.
"""
struct BnbOptions
    lb_solver::AbstractBoundingSolver
    ub_solver::AbstractBoundingSolver
    exploration::ExplorationStrategy
    branching::BranchingStrategy
    maxtime::Float64
    maxnode::Int
    tolgap::Float64
    tolint::Float64
    tolprune::Float64
    dualpruning::Bool
    l0screening::Bool
    l1screening::Bool
    verbosity::Bool
    showevery::Int
    keeptrace::Bool
    function BnbOptions(;
        lb_solver::AbstractBoundingSolver   = CDAS(),
        ub_solver::AbstractBoundingSolver   = CDAS(),
        exploration::ExplorationStrategy    = DFS,
        branching::BranchingStrategy        = LARGEST,
        maxtime::Float64                    = 60.,
        maxnode::Int                        = typemax(Int),
        tolgap::Float64                     = 1e-8,
        tolint::Float64                     = 1e-8,
        tolprune::Float64                   = 0.,
        dualpruning::Bool                   = false,
        l0screening::Bool                   = false,
        l1screening::Bool                   = false,
        verbosity::Bool                     = true,
        showevery::Int                      = 1,
        keeptrace::Bool                     = true,
    )

        @assert maxtime >= 0.
        @assert maxnode >= 0
        @assert tolgap >= 0.
        @assert tolint >= 0.
        @assert tolprune >= 0.
        @assert showevery >= 0

        return new(
            lb_solver,
            ub_solver,
            exploration,
            branching,
            maxtime,
            maxnode,
            tolgap,
            tolint,
            tolprune,
            dualpruning,
            l0screening,
            l1screening,
            verbosity,
            showevery,
            keeptrace,
        )
    end
end

mutable struct BnbNode
    parent::Union{BnbNode,Nothing}
    S0::BitArray
    S1::BitArray
    Sb::BitArray
    lb::Float64
    ub::Float64
    x::Vector
    w::Vector
    u::Vector
    x_ub::Vector
    function BnbNode(problem::Problem)
        return new(
            nothing,
            falses(problem.n),
            falses(problem.n),
            trues(problem.n),
            -Inf,
            Inf,
            zeros(problem.n),
            zeros(problem.m),
            -gradient(problem.F, problem.y, zeros(problem.m)),
            zeros(problem.n),
        )
    end
    function BnbNode(parent::BnbNode, j::Int, jval::Int, prob::Problem)
        child = new(
            parent,
            copy(parent.S0),
            copy(parent.S1),
            copy(parent.Sb),
            copy(parent.lb),
            copy(parent.ub),
            copy(parent.x),
            copy(parent.w),
            copy(parent.u),
            copy(parent.x_ub),
        )
        fixto!(child, j, jval, prob)
        return child
    end
end

"""
    BnbTrace

Trace of a [`BnbSolver`](@ref).
"""
Base.@kwdef mutable struct BnbTrace 
    ub::Vector              = Vector()
    lb::Vector              = Vector()
    node_count::Vector{Int} = Vector{Int}()
    queue_size::Vector{Int} = Vector{Int}()
    timer::Vector           = Vector()
    supp_pruned::Vector     = Vector()
    node_lb::Vector         = Vector()
    node_ub::Vector         = Vector()
    card_Sb::Vector{Int}    = Vector{Int}()
    card_S1::Vector{Int}    = Vector{Int}()
    card_S0::Vector{Int}    = Vector{Int}()
end

"""
    BnbSolver <: AbstractSolver

Branch-and-Bound algorithm to solve a [`Problem`](@ref). Additional keyword 
arguments in `kwargs` are passed to a [`BnbOptions`](@ref) instance.
"""
mutable struct BnbSolver <: AbstractSolver
    status::MOI.TerminationStatusCode
    ub::Float64
    lb::Float64
    x::Vector
    queue::Vector{BnbNode}
    node_count::Int
    supp_pruned::Float64
    start_time::Float64
    options::BnbOptions
    trace::BnbTrace
    function BnbSolver(; kwargs...)
        return new(
            MOI.OPTIMIZE_NOT_CALLED,
            Inf,
            -Inf,
            zeros(0),
            Vector{BnbNode}(),
            0,
            0.,
            Dates.time(),
            BnbOptions(; kwargs...),
            BnbTrace(),
        )
    end
end

"""
    BnbResult

Result of a [`BnbSolver`](@ref).
"""
struct BnbResult <: AbstractResult
    termination_status::MOI.TerminationStatusCode
    solve_time::Float64
    node_count::Int
    objective_value::Float64
    relative_gap::Float64
    x::Vector
    trace::BnbTrace
    function BnbResult(solver::BnbSolver, trace::BnbTrace)
        return new(
            solver.status,
            elapsed_time(solver),
            solver.node_count,
            solver.ub,
            gap(solver),
            solver.x,
            trace,
        )
    end
end

Base.show(io::IO, solver::BnbSolver) = print(io, "Bnb solver")

function initialize!(
    solver::BnbSolver, 
    problem::Problem, 
    x0::Union{Vector,Nothing}, 
    S0::Vector{Int}, 
    S1::Vector{Int},
    )
    solver.status = OPTIMIZE_NOT_CALLED
    solver.ub = isa(x0, Nothing) ? Inf : objective(problem, x0)
    solver.lb = -Inf
    solver.x = isa(x0, Nothing) ? zeros(problem.m) : x0
    root = BnbNode(problem)
    for i in S0
        fixto!(root, i, 0, problem)
    end
    for i in S1
        fixto!(root, i, 1, problem)
    end
    push!(solver.queue, root)
    solver.node_count = 0
    solver.supp_pruned = 0.
    solver.start_time = Dates.time()
    solver.trace = BnbTrace()
    return nothing
end

function header()
    str = (
        "  Nodes" *
        "   Time" * 
        "   Lower" *
        "   Upper" *
        " Rel gap" *
        " Abs gap" *
        "  %S0" *
        "  %S1" *
        "  %Sb"
    )
    return str
end

function display_head()
    str = header()
    println(repeat("-", length(str)))
    println(str)
    println(repeat("-", length(str)))
    return nothing
end

function display_trace(solver::BnbSolver, node)
    @printf "%7d" solver.node_count
    @printf " %6.2f" elapsed_time(solver)
    @printf " %7.2f" solver.lb
    @printf " %7.2f" solver.ub
    @printf "  %5.2f%%" 100 * gap(solver)
    @printf "   %.0e" abs(solver.ub - solver.lb)
    @printf " %3d%%" 100 * sum(node.S0) / length(node.S0)
    @printf " %3d%%" 100 * sum(node.S1) / length(node.S1)
    @printf " %3d%%" 100 * sum(node.Sb) / length(node.Sb)
    println()
end

function display_tail()
    str = header()
    println(repeat("-", length(str)))
end

depth(node::BnbNode) = sum(node.S0 .| node.S1)
elapsed_time(solver::BnbSolver) = Dates.time() - solver.start_time
gap(solver::BnbSolver) = abs(solver.ub - solver.lb) / (abs(solver.ub) + 1e-10)
gap(node::BnbNode) = abs(node.ub - node.lb) / (abs(node.ub) + 1e-10)
is_terminated(solver::BnbSolver) = (solver.status != OPTIMIZE_NOT_CALLED)

function update_status!(solver::BnbSolver, options::BnbOptions)
    if elapsed_time(solver) >= options.maxtime
        solver.status = MOI.TIME_LIMIT
    elseif solver.node_count >= options.maxnode
        solver.status = MOI.ITERATION_LIMIT
    elseif gap(solver) <= options.tolgap
        solver.status = MOI.OPTIMAL
    elseif isempty(solver.queue)
        solver.status = MOI.OPTIMAL
    end
    return (solver.status != MOI.OPTIMIZE_NOT_CALLED)
end

function next_node!(solver::BnbSolver, options::BnbOptions)
    if options.exploration == DFS
        node = pop!(solver.queue)
        solver.node_count += 1
    elseif options.exploration == BFS
        node = popfirst!(solver.queue)
        solver.node_count += 1
    else
        error("Exploration strategy '$(options.exploration) not implemented yet")
    end
    return node
end

function prune!(solver::BnbSolver, node::BnbNode, options::BnbOptions)
    pruning_test = (node.lb > solver.ub + options.tolprune)
    perfect_test = (options.tolprune <= gap(node) < options.tolgap)
    prune = (pruning_test | perfect_test)
    if prune
        solver.supp_pruned += 2. ^ (-depth(node))
    end
    return prune
end

function branch!(prob::Problem, solver::BnbSolver, node::BnbNode, options::BnbOptions)
    !any(node.Sb) && return nothing
    if options.branching == LARGEST
        jSb = argmax(abs.(node.x[node.Sb]))
    else
        error("Unsupported branching strategy $(options.branching)")
    end
    j = (1:prob.n)[node.Sb][jSb]
    node_j0 = BnbNode(node, j, 0, prob)
    node_j1 = BnbNode(node, j, 1, prob)
    push!(solver.queue, node_j0)
    push!(solver.queue, node_j1)
    return nothing
end

function fixto!(node::BnbNode, j::Int, jval::Int, prob::Problem)
    node.Sb[j] || error("Branching index $j is already fixed")
    node.Sb[j] = false
    if jval == 0
        node.S0[j] = true
        if node.x[j] != 0.0
            axpy!(-node.x[j], prob.A[:, j], node.w)
            copy!(node.u, -gradient(prob.F, prob.y, node.w))
            node.x[j] = 0.0
        end
    elseif jval == 1
        node.S1[j] = true
    end
    return nothing
end

function update_bounds!(solver::BnbSolver, node::BnbNode, options::BnbOptions)
    if (node.ub ≈ solver.ub) & (norm(node.x_ub, 0) < norm(solver.x, 0))
        solver.ub = copy(node.ub)
        solver.x = copy(node.x_ub)
        filter!(queue_node -> !prune!(solver, queue_node, options), solver.queue)
    elseif node.ub < solver.ub
        solver.ub = copy(node.ub)
        solver.x = copy(node.x_ub)
        filter!(queue_node -> !prune!(solver, queue_node, options), solver.queue)
    end
    if isempty(solver.queue)
        solver.lb = min(node.lb, solver.ub)
    else
        solver.lb = minimum([queue_node.lb for queue_node in solver.queue])
    end
end

function update_trace!(trace::BnbTrace, solver::BnbSolver, node::BnbNode, options::BnbOptions)
    push!(trace.ub, solver.ub)
    push!(trace.lb, solver.lb)
    push!(trace.node_count, solver.node_count)
    push!(trace.queue_size, length(solver.queue))
    push!(trace.timer, elapsed_time(solver))
    push!(trace.supp_pruned, solver.supp_pruned)
    push!(trace.node_lb, node.lb)
    push!(trace.node_ub, node.ub)
    push!(trace.card_Sb, sum(node.Sb))
    push!(trace.card_S0, sum(node.S0))
    push!(trace.card_S1, sum(node.S1))
    return nothing
end

"""
    optimize(
        solver::BnbSolver,
        problem::Problem;
        x0::Union{Vector,Nothing}=nothing,
        S0::Vector{Int}=Vector{Int}(),
        S1::Vector{Int}=Vector{Int}(),
    )

Optimize a [`Problem`](@ref) with a [`BnbSolver`](@ref). The argument `x0` is
used as a warm start. The arguments `S0` and `S1` can be used to impose zero
and non-zero constraints directly in the root node. They must match `x0`.
"""
function optimize(
    solver::BnbSolver,
    problem::Problem;
    x0::Union{Vector,Nothing}=nothing,
    S0::Vector{Int}=Vector{Int}(),
    S1::Vector{Int}=Vector{Int}(),
    )

    @assert !isa(problem.G, ZeroPenalty)
    !isa(x0, Nothing) && @assert (length(x0) == problem.n)
    !isa(x0, Nothing) && @assert all(x0[S0] .== 0.)
    !isa(x0, Nothing) && @assert all(x0[S1] .!= 0.)
    initialize!(solver, problem, x0, S0, S1)

    options = solver.options
    trace = solver.trace

    options.verbosity && display_head()
    while true
        update_status!(solver, options) 
        is_terminated(solver) && break
        node = next_node!(solver, options)
        bound!(options.lb_solver, problem, solver, node, options, LOWER)
        if !(prune!(solver, node, options))
            bound!(options.ub_solver, problem, solver, node, options, UPPER)
            branch!(problem, solver, node, options)
        end
        update_bounds!(solver, node, options)
        options.keeptrace && update_trace!(trace, solver, node, options)
        options.verbosity && solver.node_count % options.showevery == 0 && display_trace(solver, node)
    end
    options.verbosity && display_tail()

    return BnbResult(solver, trace)
end
