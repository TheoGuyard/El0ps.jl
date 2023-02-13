"""
    DirectResult <: AbstractResult

Result of a [`DirectSolver`](@ref).
"""
struct DirectResult <: AbstractResult
    termination_status::MOI.TerminationStatusCode
    solve_time::Float64
    node_count::Int
    objective_value::Float64
    relative_gap::Float64
    x::Vector
    function DirectResult(model::JuMP.Model)
        return new(
            JuMP.termination_status(model),
            JuMP.solve_time(model),
            try JuMP.node_count(model) catch; 0 end,
            has_values(model) ? JuMP.objective_value(model) : Inf,
            has_values(model) ? JuMP.relative_gap(model) : Inf,
            has_values(model) ? JuMP.value.(model[:x]) : zeros(length(model[:x])),
        )
    end
end

"""
    DirectSolver <: AbstractSolver

Formulate a [`Problem`](@ref) as a Mixed-Integer-Program and solve it using the
specified `optimizer`. 

# Arguments 

- `optimizer` : The `MOI.AbstractOptimizer` of the specified optimizer (example 
: `GLPK.Optimizer`).
- `options::Dict` : Options name and value to pass to the 
`optimizer` with the method `set_optimizer_attribute` from 
[`MathOptInterface`](https://github.com/jump-dev/MathOptInterface.jl).
"""
struct DirectSolver <: AbstractSolver
    optimizer
    options::Dict
    function DirectSolver(optimizer; 
        options::Dict=Dict()
    )
        return new(optimizer, options)
    end
end

Base.show(io::IO, solver::DirectSolver) = print(io, "Direct solver")

"""
    initialize_model(problem::Problem, solver::DirectSolver, x0::Vector)

Initialize the [`JuMP`](https://github.com/jump-dev/JuMP.jl) model 
```
min fcost + λ Gcost
st  Ax = w
    x ∈ R^n
    z ∈ {0,1}^n
    w ∈ R^m
    fcost ∈ R
    hcost ∈ R
```   
and set the initial value of `x` to `x0`. To complete the model, one has to 
construct the epigraph formulations of the functions `f` and `h` in the 
[`Problem`](@ref) using [`bind_model!`](@ref) and the scalar values `fcost` and
`hcost` that represent the epigraph value.
"""
function initialize_model(
    problem::Problem, 
    solver::DirectSolver, 
    x0::Union{Vector,Nothing},
    S0::Vector{Int},
    S1::Vector{Int},
    )

    model = Model(
        optimizer_with_attributes(
            solver.optimizer,
            solver.options...
        )
    )
    
    A = problem.A
    y = problem.y
    λ = problem.λ
    m = problem.m
    n = problem.n

    @variable(model, x[1:n])
    @variable(model, z[1:n], Bin)
    @variable(model, w[1:m])
    @constraint(model, A * x .== w)
    @variable(model, fcost)
    @variable(model, Gcost)
    @objective(model, Min, fcost + λ * Gcost)
    if !isa(x0, Nothing)
        for i in eachindex(x, x0)
            set_start_value(x[i], x0[i])
        end
    end
    for i in S0
        @constraint(model, z[i] == 0)
    end
    for i in S1
        @constraint(model, z[i] == 1)
    end

    return model
end

"""
    optimize(
        solver::DirectSolver,
        problem::Problem;
        x0::Union{Vector,Nothing}=nothing,
    )

Optimize a [`Problem`](@ref) with a [`DirectSolver`](@ref). The argument `x0` is
used as a warm start. The arguments `S0` and `S1` can be used to impose zero
and non-zero constraints directly in the root node. They must match `x0`.
"""
function optimize(
    solver::DirectSolver,
    problem::Problem;
    x0::Union{Vector,Nothing}=nothing,
    S0::Vector{Int}=Vector{Int}(),
    S1::Vector{Int}=Vector{Int}(),
    )
    !isa(x0, Nothing) && @assert (length(x0) == problem.n)
    !isa(x0, Nothing) && @assert all(x0[S0] .== 0.)
    !isa(x0, Nothing) && @assert all(x0[S1] .!= 0.)
    model = initialize_model(problem, solver, x0, S0, S1)
    bind_model!(problem.f, problem.y, model)
    bind_model!(problem.h, model)
    optimize!(model)
    return DirectResult(model)
end
