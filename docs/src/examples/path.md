## Fit a regularization path

> This examples shows how to fit a regularization [`Path`](@ref) associated to a [`Problem`](@ref).

First, import all the necessary modules.

```julia
using El0ps
```

Then, define some parameters and generate a synthetic instance.

```julia
k, m, n, ρ, s = 5, 10, 30, 0.1, 10.
x, A, y = synthetic_data_regression(k, m, n, ρ, s)
```

Instantiate a Least-squares data-fidelity function and a Big-M constraint with `M=1` as follows.

```julia
F = LeastSquares()
G = Bigm(1.)
```

Set the solver to use when solving a [`Problem`](@ref) instance at some value of `λ` in the regularization [`Path`](@ref).

```julia
solver = BnbSolver()
```

Finally, fit the regularization path.

```julia
path = fit_path(solver, F, G, A, y)
println(path)
```