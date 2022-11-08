using El0ps
using CPLEX

# Define data parameters
k = 5
m = 10
n = 30
ρ = 0.1
s = 10.

# Sample synthetic data
xtrue, A, y = synthetic_data_regression(k, m, n, ρ, s)
F = LeastSquares()
G = Bigm(1.)
λ = 0.1 * compute_λmax(F, G, A, y)

# Set a problem instance
problem = Problem(F, G, A, y, λ)

# Solve the problem using the Bnb method
solver = BnbSolver(verbosity=false)
result = optimize(solver, problem)
println(result)

# Solve the problem using a direct method
solver = DirectSolver(CPLEX.Optimizer, verbosity=false)
result = optimize(solver, problem)
println(result)
