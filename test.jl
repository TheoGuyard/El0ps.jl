using El0ps

k = 5
m = 10
n = 30
ρ = 0.1
s = 10.

xtrue, A, y = synthetic_data_regression(k, m, n, ρ, s)

F = LeastSquares()
G = Bigm(1.)
λ = 0.1 * compute_λmax(F, G, A, y)

problem = Problem(F, G, A, y, λ)

solver = El0ps.Solver(verbosity=true)

result = optimize(solver, problem)
print(result)


# setup = OrderedDict(
#     "dataset_name"  => "breast-cancer",
#     "datafit_name"  => "logistic",
#     "penalty_name"  => "bigm",
#     "solver_name"   => "bnb-dpl1l0",
#     "γmax"          => 1e-0,
#     "γmin"          => 1e-2,
#     "γnum"          => 11,
# )

# A, y = get_realworld_data(setup)

# println("Calibrating parameters...")
# F, G, λ, x0 = calibrate_l0learn(setup["datafit_name"], setup["penalty_name"], A, y)

# println("Solving the problem...")
# solver  = get_solver(setup["solver_name"])
# γstep   = (log10(setup["γmin"]) - log10(setup["γmax"])) / (setup["γnum"] - 1)
# γgrid   = 10 .^ (log10(setup["γmax"]):γstep:log10(setup["γmin"]))
# λmax    = compute_λmax(F, G, A, y)
# x0      = zeros(size(A)[2])
# results = Dict()
# for (i, γ) in enumerate(γgrid)
#     println("γ : $γ")
#     problem = Problem(F, G, A, y, γ * λmax)
#     result  = optimize(solver, problem, x0=x0)
#     println(result)
#     copy!(x0, result.x)
# end

# # λ = 0.1 * compute_λmax(F, G, A, y)

# # problem = Problem(F, G, A, y, λ)
# # println(problem)

# # result = optimize(DirectSolver(SCIP.Optimizer), problem)
# # println(result)
# # result = optimize(BnbSolver(dualpruning=true, l1screening=true, l0screening=true), problem)
# # println(result)