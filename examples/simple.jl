using El0ps

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

# Solver the problem
problem = Problem(F, G, A, y, λ)
solver = Solver(verbosity=true)
result = optimize(solver, problem)
print(result)
