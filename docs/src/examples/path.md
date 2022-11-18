## Fit a regularization path

This examples shows how to fit a regularization path of an L0-regularized problem.

First, import all the necessary modules.

```julia
using El0ps
```

Then, define some parameters and generate a synthetic instance.

```julia
k, m, n, ρ, s = 5, 10, 30, 0.1, 10.
x, A, y = synthetic_data_regression(k, m, n, ρ, s)
```

Instantiate a Least-Squares datafit and the Big-M constraint with $M=1$ as follows.

```julia
F = LeastSquares()
G = Bigm(1.)
```

Set the solver used to fit the regularization path.

```julia
solver = BnbSolver()
```

Finally, fit the regularization path.

```julia
path = fit_path(solver, F, G, A, y)
```