using El0ps

# Define data parameters
k = 2
m = 10
n = 30
ρ = 0.1
s = 10.

# Sample synthetic data
xtrue, A, y = synthetic_data_regression(k, m, n, ρ, s)
F = LeastSquares()
G = L2norm(1.)

# Set a solver
solver = BnbSolver(dualpruning=true, l0screening=true, l1screening=true)

# Fit an L0-regularization path
path = fit_path(solver, F, G, A, y, verbosity=true, compute_cv=true)