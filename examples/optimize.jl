using El0ps
using SCIP

# Define data parameters
k = 5
m = 10
n = 30
ρ = 0.1
s = 10.

# Sample synthetic data
xtrue, A, y = synthetic_data_regression(k, m, n, ρ, s)

# Set the loss and the Big-M constraint
F = LeastSquares()
G = Bigm(1.)

# Set the regularization strenght
λmax = compute_λmax(F, G, A, y) # value such that 0 is the solution
λ = 0.1 * λmax

# Instanciate a problem
problem = Problem(F, G, A, y, λ)

# Solve the problem using the Bnb method
solver = BnbSolver(
    dualpruning = true, 
    l0screening = true, 
    l1screening = true,
    verbosity   = false
)
result = optimize(solver, problem)
println(result)

# Solve the problem using a direct method with the SCIP optimizer with options.
solver = DirectSolver(
    SCIP.Optimizer, 
    options = Dict(
        "display/verblevel" => 0, 
        "limits/gap"        => 1e-4,
    )
)
result = optimize(solver, problem)
println(result)
