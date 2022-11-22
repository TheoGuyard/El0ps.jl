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
min Fcost + λ Ωcost
st  Ax = w
    x ∈ R^n
    z ∈ {0,1}^n
    w ∈ R^m
    Fcost ∈ R
    Gcost ∈ R
```   
and set the initial value of `x` to `x0`. To complete the model, one has to 
construct the epigraph formulations of the functions `F` and `G` in the 
[`Problem`](@ref) using [`bind_model!`](@ref) and the scalar values `Fcost` and
`Gcost` that represent the epigraph value.
"""
function initialize_model(problem::Problem, solver::DirectSolver, x0::Vector)

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
    @variable(model, Fcost)
    @variable(model, Ωcost)
    @objective(model, Min, Fcost + λ * Ωcost)
    for i in eachindex(x, x0)
        set_start_value(x[i], x0[i])
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
used as a warm start. 
"""
function optimize(
    solver::DirectSolver,
    problem::Problem;
    x0::Union{Vector,Nothing}=nothing,
    )
    x0 = isa(x0, Nothing) ? zeros(problem.n) : x0
    @assert length(x0) == problem.n
    model = initialize_model(problem, solver, x0)
    bind_model!(problem.F, problem.y, model)
    bind_model!(problem.G, model)
    optimize!(model)
    return DirectResult(model)
end
