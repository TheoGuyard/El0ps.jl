struct DirectResult <: AbstractResult
    termination_status::MOI.TerminationStatusCode
    solve_time::Float64
    node_count::Int
    objective_value::Float64
    relative_gap::Float64
    x::Vector{Float64}
    function DirectResult(model::JuMP.Model)
        return new(
            JuMP.termination_status(model),
            JuMP.solve_time(model),
            JuMP.node_count(model),
            has_values(model) ? JuMP.objective_value(model) : Inf,
            has_values(model) ? JuMP.relative_gap(model) : Inf,
            has_values(model) ? JuMP.value.(model[:x]) : zeros(length(model[:x])),
        )
    end
end

struct DirectSolver <: AbstractSolver
    optimizer
    options::Dict
    function DirectSolver(optimizer; options::Dict=Dict())
        return new(optimizer, options)
    end
end

Base.show(io::IO, solver::DirectSolver) = print(io, "Direct solver")

function initialize_model(problem::Problem, solver::DirectSolver, x0::Vector{Float64})

    model = Model(solver.optimizer)
    
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

function optimize(
    solver::DirectSolver,
    problem::Problem;
    x0::Union{Vector{Float64},Nothing}=nothing,
    )
    x0 = isa(x0, Nothing) ? zeros(problem.n) : x0
    @assert length(x0) == problem.n
    model = initialize_model(problem, solver, x0)
    for (option_name, option_value) in solver.options
        set_optimizer_attribute(model, String(option_name), option_value)
    end
    bind_model!(problem.F, problem.y, model)
    bind_model!(problem.G, model)
    optimize!(model)
    return DirectResult(model)
end